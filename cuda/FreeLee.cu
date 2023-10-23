/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University
  Copyright Equnior ASA

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
#include <math.h>
#include <stdio.h>
#include <cuda_profiler_api.h>

#define STOKES

#define NBLOCKS 1024
#define NTHREADS 256

__global__ void dvc_ScaLBL_D3Q19_FreeLeeModel_TwoFluid_Init(double *gqbar, double *mu_phi, double *ColorGrad, double Fx, double Fy, double Fz, int Np)
{
	int n;
    double p = 1.0;//NOTE: take initial pressure p=1.0
    double chem;
    double cg_x,cg_y,cg_z;

	//for (n=0; n<Np; n++){
	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
	    //........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;

		if ( n<Np ){
			chem = mu_phi[n];//chemical potential
			cg_x = ColorGrad[0*Np+n];
			cg_y = ColorGrad[1*Np+n];
			cg_z = ColorGrad[2*Np+n];

			gqbar[0*Np+n]  = 0.3333333333333333*p;
			gqbar[1*Np+n]  = 0.055555555555555555*(p - 0.5*(chem*cg_x+Fx));		//double(100*n)+1.f;
			gqbar[2*Np+n]  = 0.055555555555555555*(p - 0.5*(-chem*cg_x-Fx));	//double(100*n)+2.f;
			gqbar[3*Np+n]  = 0.055555555555555555*(p - 0.5*(chem*cg_y+Fy));	//double(100*n)+3.f;
			gqbar[4*Np+n]  = 0.055555555555555555*(p - 0.5*(-chem*cg_y-Fy));	//double(100*n)+4.f;
			gqbar[5*Np+n]  = 0.055555555555555555*(p - 0.5*(chem*cg_z+Fz));	//double(100*n)+5.f;
			gqbar[6*Np+n]  = 0.055555555555555555*(p - 0.5*(-chem*cg_z-Fz));	//double(100*n)+6.f;

			gqbar[7*Np+n]  = 0.0277777777777778*(p-0.5*(chem*(cg_x+cg_y)+Fx+Fy));   //double(100*n)+7.f;
			gqbar[8*Np+n]  = 0.0277777777777778*(p-0.5*(chem*(-cg_x-cg_y)-Fx-Fy));   //double(100*n)+8.f;
			gqbar[9*Np+n]  = 0.0277777777777778*(p-0.5*(chem*(cg_x-cg_y)+Fx-Fy));   //double(100*n)+9.f;
			gqbar[10*Np+n] = 0.0277777777777778*(p-0.5*(chem*(-cg_x+cg_y)-Fx+Fy));  //double(100*n)+10.f;

			gqbar[11*Np+n] = 0.0277777777777778*(p-0.5*(chem*(cg_x+cg_z)+Fx+Fz));  //double(100*n)+11.f;
			gqbar[12*Np+n] = 0.0277777777777778*(p-0.5*(chem*(-cg_x-cg_z)-Fx-Fz));  //double(100*n)+12.f;
			gqbar[13*Np+n] = 0.0277777777777778*(p-0.5*(chem*(cg_x-cg_z)+Fx-Fz));  //double(100*n)+13.f;
			gqbar[14*Np+n] = 0.0277777777777778*(p-0.5*(chem*(-cg_x+cg_z)-Fx+Fz));  //double(100*n)+14.f;

			gqbar[15*Np+n] = 0.0277777777777778*(p-0.5*(chem*(cg_y+cg_z)+Fy+Fz)); ;  //double(100*n)+15.f;
			gqbar[16*Np+n] = 0.0277777777777778*(p-0.5*(chem*(-cg_y-cg_z)-Fy-Fz));;  //double(100*n)+16.f;
			gqbar[17*Np+n] = 0.0277777777777778*(p-0.5*(chem*(cg_y-cg_z)+Fy-Fz)); ;  //double(100*n)+17.f;
			gqbar[18*Np+n] = 0.0277777777777778*(p-0.5*(chem*(-cg_y+cg_z)-Fy+Fz));;  //double(100*n)+18.f;
		}
	}
}

__global__ void dvc_ScaLBL_D3Q19_FreeLeeModel_SingleFluid_Init(double *gqbar, double Fx, double Fy, double Fz, int Np)
{
	int n;
    double p = 1.0;//NOTE: take initial pressure p=1.0

    //	for (n=0; n<Np; n++){
    int S = Np/NBLOCKS/NTHREADS + 1;
    for (int s=0; s<S; s++){
    	//........Get 1-D index for this thread....................
    	n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;

    	if ( n<Np ){        
    		gqbar[0*Np+n]  = 0.3333333333333333*p;
    		gqbar[1*Np+n]  = 0.055555555555555555*(p - 0.5*(Fx));		//double(100*n)+1.f;
    		gqbar[2*Np+n]  = 0.055555555555555555*(p - 0.5*(-Fx));	//double(100*n)+2.f;
    		gqbar[3*Np+n]  = 0.055555555555555555*(p - 0.5*(Fy));	//double(100*n)+3.f;
    		gqbar[4*Np+n]  = 0.055555555555555555*(p - 0.5*(-Fy));	//double(100*n)+4.f;
    		gqbar[5*Np+n]  = 0.055555555555555555*(p - 0.5*(Fz));	//double(100*n)+5.f;
    		gqbar[6*Np+n]  = 0.055555555555555555*(p - 0.5*(-Fz));	//double(100*n)+6.f;

    		gqbar[7*Np+n]  = 0.0277777777777778*(p-0.5*(Fx+Fy));   //double(100*n)+7.f;
    		gqbar[8*Np+n]  = 0.0277777777777778*(p-0.5*(-Fx-Fy));   //double(100*n)+8.f;
    		gqbar[9*Np+n]  = 0.0277777777777778*(p-0.5*(Fx-Fy));   //double(100*n)+9.f;
    		gqbar[10*Np+n] = 0.0277777777777778*(p-0.5*(-Fx+Fy));  //double(100*n)+10.f;

    		gqbar[11*Np+n] = 0.0277777777777778*(p-0.5*(Fx+Fz));  //double(100*n)+11.f;
    		gqbar[12*Np+n] = 0.0277777777777778*(p-0.5*(-Fx-Fz));  //double(100*n)+12.f;
    		gqbar[13*Np+n] = 0.0277777777777778*(p-0.5*(Fx-Fz));  //double(100*n)+13.f;
    		gqbar[14*Np+n] = 0.0277777777777778*(p-0.5*(-Fx+Fz));  //double(100*n)+14.f;

    		gqbar[15*Np+n] = 0.0277777777777778*(p-0.5*(Fy+Fz)); ;  //double(100*n)+15.f;
    		gqbar[16*Np+n] = 0.0277777777777778*(p-0.5*(-Fy-Fz));;  //double(100*n)+16.f;
    		gqbar[17*Np+n] = 0.0277777777777778*(p-0.5*(Fy-Fz)); ;  //double(100*n)+17.f;
    		gqbar[18*Np+n] = 0.0277777777777778*(p-0.5*(-Fy+Fz));;  //double(100*n)+18.f;
    	}
    }
}

__global__ void dvc_ScaLBL_FreeLeeModel_PhaseField_Init(int *Map, double *Phi, double *Den, double *hq, double *ColorGrad, 
                                                    double rhoA, double rhoB, double tauM, double W, int start, int finish, int Np){
	int idx,n;
	double phi;
    double nx,ny,nz,cg_mag;
    double theta;//anti-diffusion term
    double cs2_inv = 4.5;//inverse of speed of sound for D3Q7
    double M = 1.0/cs2_inv*(tauM-0.5);//diffusivity (or mobility)

	for (idx=start; idx<finish; idx++){

		n = Map[idx];
		phi = Phi[n];
		if (phi > 1.f)   phi = 1.0;
		if (phi < -1.f)  phi = -1.0;
		Den[idx] = rhoA + 0.5*(1.0-phi)*(rhoB-rhoA);

        //compute unit normal of color gradient
        nx = ColorGrad[idx+0*Np];
        ny = ColorGrad[idx+1*Np];
        nz = ColorGrad[idx+2*Np];
        cg_mag = sqrt(nx*nx+ny*ny+nz*nz);
		double ColorMag_temp = cg_mag;
		if (cg_mag==0.0) ColorMag_temp=1.0;
		nx = nx/ColorMag_temp;
		ny = ny/ColorMag_temp;
		nz = nz/ColorMag_temp;		

        theta = M*cs2_inv*2.0*(1-phi*phi)/W;
		hq[0*Np+idx]=0.3333333333333333*(phi);
		hq[1*Np+idx]=0.1111111111111111*(phi+theta*nx);
		hq[2*Np+idx]=0.1111111111111111*(phi-theta*nx);
		hq[3*Np+idx]=0.1111111111111111*(phi+theta*ny);
		hq[4*Np+idx]=0.1111111111111111*(phi-theta*ny);
		hq[5*Np+idx]=0.1111111111111111*(phi+theta*nz);
		hq[6*Np+idx]=0.1111111111111111*(phi-theta*nz);
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAodd_FreeLeeModel_PhaseField(int *neighborList, int *Map, double *hq, double *Den, double *Phi, 
                                                          double rhoA, double rhoB, int start, int finish, int Np){

	int idx,n,nread;
	double fq,phi;
	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

		if ( n<finish ){

			// q=0
			fq = hq[n];
			phi = fq;

			// q=1
			nread = neighborList[n]; 
			fq = hq[nread];
			phi += fq;

			// q=2
			nread = neighborList[n+Np]; 
			fq = hq[nread];  
			phi += fq;

			// q=3
			nread = neighborList[n+2*Np]; 
			fq = hq[nread];
			phi += fq;

			// q = 4
			nread = neighborList[n+3*Np]; 
			fq = hq[nread];
			phi += fq;

			// q=5
			nread = neighborList[n+4*Np];
			fq = hq[nread];
			phi += fq;

			// q = 6
			nread = neighborList[n+5*Np];
			fq = hq[nread];
			phi += fq;

			// save the number densities
			Den[n] = rhoA + 0.5*(1.0-phi)*(rhoB-rhoA);

			// save the phase indicator field
			idx = Map[n];
			Phi[idx] = phi; 
		}
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAeven_FreeLeeModel_PhaseField( int *Map, double *hq, double *Den, double *Phi,
		double rhoA, double rhoB, int start, int finish, int Np){
	
	int idx,n;
	double fq, phi;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

		if ( n<finish ){

            // q=0
            fq = hq[n];
            phi = fq;
            
            // q=1
            fq = hq[2*Np+n];
            phi += fq;

            // f2 = hq[10*Np+n];
            fq = hq[1*Np+n];
            phi += fq;

            // q=3
            fq = hq[4*Np+n];
            phi += fq;

            // q = 4
            fq = hq[3*Np+n];
            phi += fq;

            // q=5
            fq = hq[6*Np+n];
            phi += fq;

            // q = 6
            fq = hq[5*Np+n];
            phi += fq;

            // save the number densities
            Den[n] = rhoA + 0.5*(1.0-phi)*(rhoB-rhoA);

            // save the phase indicator field
            idx = Map[n];
            Phi[idx] = phi; 	
		}
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAodd_FreeLee_PhaseField(int *neighborList, int *Map, double *hq, double *Den, double *Phi, double *ColorGrad, double *Vel,
		double rhoA, double rhoB, double tauM, double W, int start, int finish, int Np){

	int n,idx,nr1,nr2,nr3,nr4,nr5,nr6;
	double h0,h1,h2,h3,h4,h5,h6;
	double nx,ny,nz,C;
	double ux,uy,uz;
	double phi, theta;
	double M = 2.0/9.0*(tauM-0.5);//diffusivity (or mobility) for the phase field D3Q7

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

		if ( n<finish ){		

			/* load phase indicator field */
			idx = Map[n];
			phi = Phi[idx];
            theta = 4.5*M*2.0*(1-phi*phi)/W;

			/* velocity */
			ux = Vel[0*Np+n];
			uy = Vel[1*Np+n];
			uz = Vel[2*Np+n];

	        /*color gradient */
			nx = ColorGrad[0*Np+n];
			ny = ColorGrad[1*Np+n];
			nz = ColorGrad[2*Np+n];
			
			//Normalize the Color Gradient
			C = sqrt(nx*nx+ny*ny+nz*nz);
			double ColorMag = C;
			if (C < 1.0e-12) ColorMag=1.0;
			nx = nx/ColorMag;
			ny = ny/ColorMag;
			nz = nz/ColorMag;		

			// q=1
			nr1 = neighborList[n]; 
			nr2 = neighborList[n+Np]; 
			nr3 = neighborList[n+2*Np]; 
			nr4 = neighborList[n+3*Np];
			nr5 = neighborList[n+4*Np];
			nr6 = neighborList[n+5*Np];
			
			//q=0
			h0 = hq[n];
			//q=1
			h1 = hq[nr1]; 
			//q=2
			h2 = hq[nr2];  
			//q=3
			h3 = hq[nr3];
			//q=4
			h4 = hq[nr4];
			//q=5
			h5 = hq[nr5];
			//q=6
			h6 = hq[nr6];

	        //-------------------------------- BGK collison for phase field ---------------------------------//
			h0 -= (h0 - 0.3333333333333333*phi)/tauM;
            h1 -= (h1 - 0.1111111111111111*nx*theta - phi*(0.1111111111111111 + 0.5*ux))/tauM;
            h2 -= (h2 + 0.1111111111111111*nx*theta - phi*(0.1111111111111111 - 0.5*ux))/tauM;
            h3 -= (h3 - 0.1111111111111111*ny*theta - phi*(0.1111111111111111 + 0.5*uy))/tauM;
            h4 -= (h4 + 0.1111111111111111*ny*theta - phi*(0.1111111111111111 - 0.5*uy))/tauM;
            h5 -= (h5 - 0.1111111111111111*nz*theta - phi*(0.1111111111111111 + 0.5*uz))/tauM;
            h6 -= (h6 + 0.1111111111111111*nz*theta - phi*(0.1111111111111111 - 0.5*uz))/tauM; 
			//........................................................................
			
			/*Update the distributions */
			// q = 0
			hq[n] = h0;
			hq[nr2] = h1;
			hq[nr1] = h2;
			hq[nr4] = h3;
			hq[nr3] = h4;
			hq[nr6] = h5;
			hq[nr5] = h6;
			//........................................................................
		}
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAeven_FreeLee_PhaseField( int *Map, double *hq, double *Den, double *Phi, double *ColorGrad, double *Vel,
		double rhoA, double rhoB, double tauM, double W, int start, int finish, int Np){
	
	int idx,n;
	double h0,h1,h2,h3,h4,h5,h6;
	double nx,ny,nz,C;
	double ux,uy,uz;
	double phi, theta;
    double M = 2.0/9.0*(tauM-0.5);//diffusivity (or mobility) for the phase field D3Q7

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

		if ( n<finish ){
			/* load phase indicator field */
			idx = Map[n];
			phi = Phi[idx];
	        theta = 4.5*M*2.0*(1-phi*phi)/W;

			/* velocity */
			ux = Vel[0*Np+n];
			uy = Vel[1*Np+n];
			uz = Vel[2*Np+n];

			/*color gradient */
			nx = ColorGrad[0*Np+n];
			ny = ColorGrad[1*Np+n];
			nz = ColorGrad[2*Np+n];
			//Normalize the Color Gradient
			C = sqrt(nx*nx+ny*ny+nz*nz);
			double ColorMag = C;
			if (C < 1.0e-12) ColorMag=1.0;
			nx = nx/ColorMag;
			ny = ny/ColorMag;
			nz = nz/ColorMag;

			h0 = hq[n];
			h1 = hq[2*Np+n]; 
			h2 = hq[Np+n];  
			h3 = hq[4*Np+n];
			h4 = hq[3*Np+n];
			h5 = hq[6*Np+n];
			h6 = hq[5*Np+n];

			//-------------------------------- BGK collison for phase field ---------------------------------//
			h0 -= (h0 - 0.3333333333333333*phi)/tauM;
            h1 -= (h1 - 0.1111111111111111*nx*theta - phi*(0.1111111111111111 + 0.5*ux))/tauM;
            h2 -= (h2 + 0.1111111111111111*nx*theta - phi*(0.1111111111111111 - 0.5*ux))/tauM;
            h3 -= (h3 - 0.1111111111111111*ny*theta - phi*(0.1111111111111111 + 0.5*uy))/tauM;
            h4 -= (h4 + 0.1111111111111111*ny*theta - phi*(0.1111111111111111 - 0.5*uy))/tauM;
            h5 -= (h5 - 0.1111111111111111*nz*theta - phi*(0.1111111111111111 + 0.5*uz))/tauM;
            h6 -= (h6 + 0.1111111111111111*nz*theta - phi*(0.1111111111111111 - 0.5*uz))/tauM; 
			//........................................................................

			/*Update the distributions */
			// q = 0
			hq[n] = h0;
			hq[Np+n] = h1;
			hq[2*Np+n] = h2;
			hq[3*Np+n] = h3;
			hq[4*Np+n] = h4;
			hq[5*Np+n] = h5;
			hq[6*Np+n] = h6;
			//........................................................................
			
		}
	}
}

__global__ void dvc_ScaLBL_D3Q7_ComputePhaseField(int *Map,  double *hq, double *Den, double *Phi, double rhoA, double rhoB, int start, int finish, int Np){
	int idx,n;
	double h0,h1,h2,h3,h4,h5,h6;
	double phi;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

		if ( n<finish ){
			h0 = hq[n];
			h1 = hq[1*Np+n]; 
			h2 = hq[2*Np+n];  
			h3 = hq[3*Np+n];
			h4 = hq[4*Np+n];
			h5 = hq[5*Np+n];
			h6 = hq[6*Np+n];

			phi = h0+h1+h2+h3+h4+h5+h6;

			// save the number densities
			Den[n] = rhoA + 0.5*(1.0-phi)*(rhoB-rhoA);

			// save the phase indicator field
			idx = Map[n];
			Phi[idx] = phi; 
		}
	}
}

__global__ void dvc_ScaLBL_D3Q19_AAodd_FreeLeeModel(int *neighborList, int *Map, double *dist, double *Den,	double *Phi, double *mu_phi, double *Vel, double *Pressure, double *ColorGrad, 
        double rhoA, double rhoB, double tauA, double tauB, double kappa, double beta, double W, double Fx, double Fy, double Fz, 
        int strideY, int strideZ, int start, int finish, int Np){

	int n,nn,nn2x,ijk;
	int nr1,nr2,nr3,nr4,nr5,nr6,nr7,nr8,nr9,nr10,nr11,nr12,nr13,nr14,nr15,nr16,nr17,nr18;
    double ux,uy,uz;//fluid velocity 
    double p;//pressure
    double chem;//chemical potential
    double phi; //phase field
    double rho0;//fluid density
	// distribution functions
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
	double m0,m3,m5,m7;
	double mm1,mm2,mm4,mm6,mm8,mm9,mm10,mm11,mm12,mm13,mm14,mm15,mm16,mm17,mm18;
	double mm3,mm5,mm7;
    double feq0,feq1,feq2,feq3,feq4,feq5,feq6,feq7,feq8,feq9,feq10,feq11,feq12,feq13,feq14,feq15,feq16,feq17,feq18;
    double nx,ny,nz;//normal color gradient
    double mgx,mgy,mgz;//mixed gradient reaching secondary neighbor

	//double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18;
    //double h0,h1,h2,h3,h4,h5,h6;//distributions for LB phase field
	double tau;//position dependent LB relaxation time for fluid
    //double C,theta;
    // double M = 2.0/9.0*(tauM-0.5);//diffusivity (or mobility) for the phase field D3Q7
    //	for (int n=start; n<finish; n++){
    int S = Np/NBLOCKS/NTHREADS + 1;
    for (int s=0; s<S; s++){
    	//........Get 1-D index for this thread....................
    	n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

    	if ( n<finish ){
    		rho0 = Den[n];//load density

    		// Get the 1D index based on regular data layout
    		ijk = Map[n];
            phi = Phi[ijk];// load phase field

    		// local relaxation time
    		tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);

    		//					COMPUTE THE COLOR GRADIENT
    		//........................................................................
    		//.................Read Phase Indicator Values............................
    		//........................................................................
    		nn = ijk-1;							// neighbor index (get convention)
    		m1 = Phi[nn];						// get neighbor for phi - 1
    		//........................................................................
    		nn = ijk+1;							// neighbor index (get convention)
    		m2 = Phi[nn];						// get neighbor for phi - 2
    		//........................................................................
    		nn = ijk-strideY;							// neighbor index (get convention)
    		m3 = Phi[nn];					// get neighbor for phi - 3
    		//........................................................................
    		nn = ijk+strideY;							// neighbor index (get convention)
    		m4 = Phi[nn];					// get neighbor for phi - 4
    		//........................................................................
    		nn = ijk-strideZ;						// neighbor index (get convention)
    		m5 = Phi[nn];					// get neighbor for phi - 5
    		//........................................................................
    		nn = ijk+strideZ;						// neighbor index (get convention)
    		m6 = Phi[nn];					// get neighbor for phi - 6
    		//........................................................................
    		nn = ijk-strideY-1;						// neighbor index (get convention)
    		m7 = Phi[nn];					// get neighbor for phi - 7
    		//........................................................................
    		nn = ijk+strideY+1;						// neighbor index (get convention)
    		m8 = Phi[nn];					// get neighbor for phi - 8
    		//........................................................................
    		nn = ijk+strideY-1;						// neighbor index (get convention)
    		m9 = Phi[nn];					// get neighbor for phi - 9
    		//........................................................................
    		nn = ijk-strideY+1;						// neighbor index (get convention)
    		m10 = Phi[nn];					// get neighbor for phi - 10
    		//........................................................................
    		nn = ijk-strideZ-1;						// neighbor index (get convention)
    		m11 = Phi[nn];					// get neighbor for phi - 11
    		//........................................................................
    		nn = ijk+strideZ+1;						// neighbor index (get convention)
    		m12 = Phi[nn];					// get neighbor for phi - 12
    		//........................................................................
    		nn = ijk+strideZ-1;						// neighbor index (get convention)
    		m13 = Phi[nn];					// get neighbor for phi - 13
    		//........................................................................
    		nn = ijk-strideZ+1;						// neighbor index (get convention)
    		m14 = Phi[nn];					// get neighbor for phi - 14
    		//........................................................................
    		nn = ijk-strideZ-strideY;					// neighbor index (get convention)
    		m15 = Phi[nn];					// get neighbor for phi - 15
    		//........................................................................
    		nn = ijk+strideZ+strideY;					// neighbor index (get convention)
    		m16 = Phi[nn];					// get neighbor for phi - 16
    		//........................................................................
    		nn = ijk+strideZ-strideY;					// neighbor index (get convention)
    		m17 = Phi[nn];					// get neighbor for phi - 17
    		//........................................................................
    		nn = ijk-strideZ+strideY;					// neighbor index (get convention)
    		m18 = Phi[nn];					// get neighbor for phi - 18

            // compute mixed difference (Eq.30, A.Fukhari et al. JCP 315(2016) 434-457)
    		//........................................................................
    		nn2x = ijk-2;							// neighbor index (get convention)
    		mm1 = Phi[nn2x];						// get neighbor for phi - 1
            mm1 = 0.25*(-mm1+5.0*m1-3.0*phi-m2);
    		//........................................................................
    		nn2x = ijk+2;							// neighbor index (get convention)
    		mm2 = Phi[nn2x];						// get neighbor for phi - 2
            mm2 = 0.25*(-mm2+5.0*m2-3.0*phi-m1);
    		//........................................................................
    		nn2x = ijk-strideY*2;							// neighbor index (get convention)
    		mm3 = Phi[nn2x];					// get neighbor for phi - 3
            mm3 = 0.25*(-mm3+5.0*m3-3.0*phi-m4);
    		//........................................................................
    		nn2x = ijk+strideY*2;							// neighbor index (get convention)
    		mm4 = Phi[nn2x];					// get neighbor for phi - 4
            mm4 = 0.25*(-mm4+5.0*m4-3.0*phi-m3);
    		//........................................................................
    		nn2x = ijk-strideZ*2;						// neighbor index (get convention)
    		mm5 = Phi[nn2x];					// get neighbor for phi - 5
            mm5 = 0.25*(-mm5+5.0*m5-3.0*phi-m6);
    		//........................................................................
    		nn2x = ijk+strideZ*2;						// neighbor index (get convention)
    		mm6 = Phi[nn2x];					// get neighbor for phi - 6
            mm6 = 0.25*(-mm6+5.0*m6-3.0*phi-m5);
    		//........................................................................
    		nn2x = ijk-strideY*2-2;						// neighbor index (get convention)
    		mm7 = Phi[nn2x];					// get neighbor for phi - 7
            mm7 = 0.25*(-mm7+5.0*m7-3.0*phi-m8);
    		//........................................................................
    		nn2x = ijk+strideY*2+2;						// neighbor index (get convention)
    		mm8 = Phi[nn2x];					// get neighbor for phi - 8
            mm8 = 0.25*(-mm8+5.0*m8-3.0*phi-m7);
    		//........................................................................
    		nn2x = ijk+strideY*2-2;						// neighbor index (get convention)
    		mm9 = Phi[nn2x];					// get neighbor for phi - 9
            mm9 = 0.25*(-mm9+5.0*m9-3.0*phi-m10);
    		//........................................................................
    		nn2x = ijk-strideY*2+2;						// neighbor index (get convention)
    		mm10 = Phi[nn2x];					// get neighbor for phi - 10
            mm10 = 0.25*(-mm10+5.0*m10-3.0*phi-m9);
    		//........................................................................
    		nn2x = ijk-strideZ*2-2;						// neighbor index (get convention)
    		mm11 = Phi[nn2x];					// get neighbor for phi - 11
            mm11 = 0.25*(-mm11+5.0*m11-3.0*phi-m12);
    		//........................................................................
    		nn2x = ijk+strideZ*2+2;						// neighbor index (get convention)
    		mm12 = Phi[nn2x];					// get neighbor for phi - 12
            mm12 = 0.25*(-mm12+5.0*m12-3.0*phi-m11);
    		//........................................................................
    		nn2x = ijk+strideZ*2-2;						// neighbor index (get convention)
    		mm13 = Phi[nn2x];					// get neighbor for phi - 13
            mm13 = 0.25*(-mm13+5.0*m13-3.0*phi-m14);
    		//........................................................................
    		nn2x = ijk-strideZ*2+2;						// neighbor index (get convention)
    		mm14 = Phi[nn2x];					// get neighbor for phi - 14
            mm14 = 0.25*(-mm14+5.0*m14-3.0*phi-m13);
    		//........................................................................
    		nn2x = ijk-strideZ*2-strideY*2;					// neighbor index (get convention)
    		mm15 = Phi[nn2x];					// get neighbor for phi - 15
            mm15 = 0.25*(-mm15+5.0*m15-3.0*phi-m16);
    		//........................................................................
    		nn2x = ijk+strideZ*2+strideY*2;					// neighbor index (get convention)
    		mm16 = Phi[nn2x];					// get neighbor for phi - 16
            mm16 = 0.25*(-mm16+5.0*m16-3.0*phi-m15);
    		//........................................................................
    		nn2x = ijk+strideZ*2-strideY*2;					// neighbor index (get convention)
    		mm17 = Phi[nn2x];					// get neighbor for phi - 17
            mm17 = 0.25*(-mm17+5.0*m17-3.0*phi-m18);
    		//........................................................................
    		nn2x = ijk-strideZ*2+strideY*2;					// neighbor index (get convention)
    		mm18 = Phi[nn2x];					// get neighbor for phi - 18
            mm18 = 0.25*(-mm18+5.0*m18-3.0*phi-m17);


    		//............Compute the Color Gradient...................................
    		nx = -3.0*1.0/18.0*(m1-m2+0.5*(m7-m8+m9-m10+m11-m12+m13-m14));
    		ny = -3.0*1.0/18.0*(m3-m4+0.5*(m7-m8-m9+m10+m15-m16+m17-m18));
    		nz = -3.0*1.0/18.0*(m5-m6+0.5*(m11-m12-m13+m14+m15-m16-m17+m18));
    		//............Compute the Chemical Potential...............................
            chem = 2.0*3.0/18.0*(m1+m2+m3+m4+m5+m6-6*phi+0.5*(m7+m8+m9+m10+m11+m12+m13+m14+m15+m16+m17+m18-12*phi));//intermediate var, i.e. the laplacian
            chem = 4.0*beta*phi*(phi+1.0)*(phi-1.0)-kappa*chem;
    		//............Compute the Mixed Gradient...................................
    		mgx = -3.0*1.0/18.0*(mm1-mm2+0.5*(mm7-mm8+mm9-mm10+mm11-mm12+mm13-mm14));
    		mgy = -3.0*1.0/18.0*(mm3-mm4+0.5*(mm7-mm8-mm9+mm10+mm15-mm16+mm17-mm18));
    		mgz = -3.0*1.0/18.0*(mm5-mm6+0.5*(mm11-mm12-mm13+mm14+mm15-mm16-mm17+mm18));
    		
    		// q=0
    		m0 = dist[n];
    		// q=1
    		nr1 = neighborList[n]; // neighbor 2 ( > 10Np => odd part of dist)
    		m1 = dist[nr1]; // reading the f1 data into register fq

    		nr2 = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
    		m2 = dist[nr2];  // reading the f2 data into register fq

    		// q=3
    		nr3 = neighborList[n+2*Np]; // neighbor 4
    		m3 = dist[nr3];

    		// q = 4
    		nr4 = neighborList[n+3*Np]; // neighbor 3
    		m4 = dist[nr4];

    		// q=5
    		nr5 = neighborList[n+4*Np];
    		m5 = dist[nr5];

    		// q = 6
    		nr6 = neighborList[n+5*Np];
    		m6 = dist[nr6];
    		
    		// q=7
    		nr7 = neighborList[n+6*Np];
    		m7 = dist[nr7];

    		// q = 8
    		nr8 = neighborList[n+7*Np];
    		m8 = dist[nr8];

    		// q=9
    		nr9 = neighborList[n+8*Np];
    		m9 = dist[nr9];

    		// q = 10
    		nr10 = neighborList[n+9*Np];
    		m10 = dist[nr10];

    		// q=11
    		nr11 = neighborList[n+10*Np];
    		m11 = dist[nr11];

    		// q=12
    		nr12 = neighborList[n+11*Np];
    		m12 = dist[nr12];

    		// q=13
    		nr13 = neighborList[n+12*Np];
    		m13 = dist[nr13];

    		// q=14
    		nr14 = neighborList[n+13*Np];
    		m14 = dist[nr14];

    		// q=15
    		nr15 = neighborList[n+14*Np];
    		m15 = dist[nr15];

    		// q=16
    		nr16 = neighborList[n+15*Np];
    		m16 = dist[nr16];

    		// q=17
    		nr17 = neighborList[n+16*Np];
    		m17 = dist[nr17];

    		// q=18
    		nr18 = neighborList[n+17*Np];
    		m18 = dist[nr18];

            //compute fluid velocity
            ux = 3.0/rho0*(m1-m2+m7-m8+m9-m10+m11-m12+m13-m14+0.5*(chem*nx+Fx)/3.0);
            uy = 3.0/rho0*(m3-m4+m7-m8-m9+m10+m15-m16+m17-m18+0.5*(chem*ny+Fy)/3.0);
            uz = 3.0/rho0*(m5-m6+m11-m12-m13+m14+m15-m16-m17+m18+0.5*(chem*nz+Fz)/3.0);
            //compute pressure
            p = (m0+m2+m1+m4+m3+m6+m5+m8+m7+m10+m9+m12+m11+m14+m13+m16+m15+m18+m17)
                      +0.5*(rhoA-rhoB)/2.0/3.0*(ux*nx+uy*ny+uz*nz);

            //compute equilibrium distributions
                feq0 = 0.3333333333333333*p - 0.25*(Fx*ux + Fy*uy + Fz*uz)*(-0.6666666666666666 + ux*ux + uy*uy + uz*uz) - 
         0.16666666666666666*rho0*(ux*ux + uy*uy + uz*uz) - 0.5*(-(nx*ux) - ny*uy - nz*uz)*
          (-0.08333333333333333*(rhoA - rhoB)*(ux*ux + uy*uy + uz*uz) + chem*(0.3333333333333333 - 0.5*(ux*ux + uy*uy + uz*uz)));
                feq1 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-ux*ux + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*uz)) - 
         0.125*(Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 
           0.3333333333333333*(-2.*ux + ux*ux + uy*uy + uz*uz)) - 0.0625*(nx - nx*ux - ny*uy - nz*uz)*
          (2*chem*ux*ux - 0.3333333333333333*((-rhoA + rhoB)*ux*ux + 2*chem*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
           0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*uz)));
                feq2 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-ux*ux + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*uz)) - 
         0.125*(Fx + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 
           0.3333333333333333*(2.*ux + ux*ux + uy*uy + uz*uz)) - 0.0625*(nx + nx*ux + ny*uy + nz*uz)*
          (-2.*chem*ux*ux + 0.1111111111111111*(-4.*chem + rhoB*(-2.*ux - 1.*ux*ux - 1.*uy*uy - 1.*uz*uz) + 
             rhoA*(2.*ux + ux*ux + uy*uy + uz*uz)) + 0.3333333333333333*((-1.*rhoA + rhoB)*ux*ux + 
             chem*(4.*ux + 2.*ux*ux + 2.*uy*uy + 2.*uz*uz)));
                feq3 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uy*uy + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*uz)) - 
         0.125*(Fx*ux + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy + 
           0.3333333333333333*(ux*ux - 2.*uy + uy*uy + uz*uz)) - 0.0625*(ny - nx*ux - ny*uy - nz*uz)*
          (2*chem*uy*uy - 0.3333333333333333*((-rhoA + rhoB)*uy*uy + 2*chem*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
           0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*uz)));
                feq4 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uy*uy + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*uz)) - 
         0.125*(Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy + 
           0.3333333333333333*(ux*ux + 2.*uy + uy*uy + uz*uz)) - 0.0625*(ny + nx*ux + ny*uy + nz*uz)*
          (-2.*chem*uy*uy + 0.1111111111111111*(-4.*chem + rhoB*(-1.*ux*ux - 2.*uy - 1.*uy*uy - 1.*uz*uz) + 
             rhoA*(ux*ux + 2.*uy + uy*uy + uz*uz)) + 0.3333333333333333*((-1.*rhoA + rhoB)*uy*uy + 
             chem*(2.*ux*ux + 4.*uy + 2.*uy*uy + 2.*uz*uz))); 
                feq5 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uz*uz + 0.3333333333333333*(ux*ux + uy*uy + (-2 + uz)*uz)) - 
         0.125*(Fx*ux + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*uz*uz + 
           0.3333333333333333*(ux*ux + uy*uy + (-2. + uz)*uz)) - 0.0625*(nx*ux + ny*uy + nz*(-1. + uz))*
          (-2.*chem*uz*uz + 0.1111111111111111*(-4.*chem + rhoB*(-1.*ux*ux - 1.*uy*uy + (2. - 1.*uz)*uz) + 
             rhoA*(ux*ux + uy*uy + (-2. + uz)*uz)) + 0.3333333333333333*((-1.*rhoA + rhoB)*uz*uz + 
             chem*(2.*ux*ux + 2.*uy*uy + uz*(-4. + 2.*uz)))); 
                feq6 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uz*uz + 0.3333333333333333*(ux*ux + uy*uy + uz*(2 + uz))) - 
         0.125*(Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*uz*uz + 
           0.3333333333333333*(ux*ux + uy*uy + uz*(2. + uz))) - 0.0625*(nz + nx*ux + ny*uy + nz*uz)*
          (-2.*chem*uz*uz + 0.1111111111111111*(-4.*chem + rhoB*(-1.*ux*ux - 1.*uy*uy + (-2. - 1.*uz)*uz) + 
             rhoA*(ux*ux + uy*uy + uz*(2. + uz))) + 0.3333333333333333*((-1.*rhoA + rhoB)*uz*uz + 
             chem*(2.*ux*ux + 2.*uy*uy + uz*(4. + 2.*uz)))); 
                feq7 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux + uy)*(ux + uy) + 0.3333333333333333*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) - 
         0.0625*(Fx*(-1. + ux) + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uy - 1.*uy*uy + 
           0.3333333333333333*(-2.*ux + ux*ux - 2.*uy + uy*uy + uz*uz)) - 0.03125*(nx + ny - nx*ux - ny*uy - nz*uz)*
          (2*chem*(ux + uy)*(ux + uy) + 0.3333333333333333*((rhoA - rhoB)*(ux + uy)*(ux + uy) - 2*chem*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
           0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz))); 
                feq8 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux + uy)*(ux + uy) + 0.3333333333333333*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) - 
         0.0625*(Fx + Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uy - 1.*uy*uy + 
           0.3333333333333333*(2.*ux + ux*ux + 2.*uy + uy*uy + uz*uz)) - 0.03125*(-(nx*(1 + ux)) - ny*(1 + uy) - nz*uz)*
          (2*chem*(ux + uy)*(ux + uy) - 0.3333333333333333*(-((rhoA - rhoB)*(ux + uy)*(ux + uy)) + 
             2*chem*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz))); 
                feq9 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux - uy)*(ux - uy) + 0.3333333333333333*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) - 
         0.0625*(Fy + Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uy - 1.*uy*uy + 
           0.3333333333333333*(-2.*ux + ux*ux + 2.*uy + uy*uy + uz*uz)) - 0.03125*(nx - nx*ux - ny*(1 + uy) - nz*uz)*
          (2*chem*(ux - uy)*(ux - uy) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uy)*(ux - uy)) + 
             2*chem*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz))); 
                feq10 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux - uy)*(ux - uy) + 0.3333333333333333*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) - 
         0.0625*(Fx*(1 + ux) + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uy - 1.*uy*uy + 
           0.3333333333333333*(2.*ux + ux*ux - 2.*uy + uy*uy + uz*uz)) - 0.03125*(ny - nx*(1 + ux) - ny*uy - nz*uz)*
          (2*chem*(ux - uy)*(ux - uy) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uy)*(ux - uy)) + 
             2*chem*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz))); 
                feq11 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux + uz)*(ux + uz) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) - 
         0.0625*(Fx*(-1. + ux) + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uz - 1.*uz*uz + 
           0.3333333333333333*(-2.*ux + ux*ux + uy*uy + (-2. + uz)*uz)) - 0.03125*(nx + nz - nx*ux - ny*uy - nz*uz)*
          (2*chem*(ux + uz)*(ux + uz) + 0.3333333333333333*((rhoA - rhoB)*(ux + uz)*(ux + uz) - 2*chem*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
           0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz))); 
                feq12 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux + uz)*(ux + uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) - 
         0.0625*(Fx + Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uz - 1.*uz*uz + 
           0.3333333333333333*(2.*ux + ux*ux + uy*uy + uz*(2. + uz))) - 0.03125*(-(nx*(1 + ux)) - ny*uy - nz*(1 + uz))*
          (2*chem*(ux + uz)*(ux + uz) - 0.3333333333333333*(-((rhoA - rhoB)*(ux + uz)*(ux + uz)) + 
             2*chem*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))); 
                feq13 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux - uz)*(ux - uz) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) - 
         0.0625*(Fz + Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uz - 1.*uz*uz + 
           0.3333333333333333*(-2.*ux + ux*ux + uy*uy + uz*(2. + uz))) - 0.03125*(nx - nx*ux - ny*uy - nz*(1 + uz))*
          (2*chem*(ux - uz)*(ux - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uz)*(ux - uz)) + 
             2*chem*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*(2 + uz)))); 
                feq14 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux - uz)*(ux - uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) - 
         0.0625*(Fx*(1 + ux) + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uz - 1.*uz*uz + 
           0.3333333333333333*(2.*ux + ux*ux + uy*uy + (-2. + uz)*uz)) - 0.03125*(nz - nx*(1 + ux) - ny*uy - nz*uz)*
          (2*chem*(ux - uz)*(ux - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uz)*(ux - uz)) + 
             2*chem*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz))); 
                feq15 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) - 
         0.0625*(Fx*ux + Fy*(-1. + uy) + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*uy*uy - 2.*uy*uz - 1.*uz*uz + 
           0.3333333333333333*(ux*ux - 2.*uy + uy*uy + (-2. + uz)*uz)) - 0.03125*(ny + nz - nx*ux - ny*uy - nz*uz)*
          (2*chem*(uy + uz)*(uy + uz) + 0.3333333333333333*((rhoA - rhoB)*(uy + uz)*(uy + uz) - 2*chem*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) + 
           0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))); 
                feq16 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) - 
         0.0625*(Fy + Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy - 2.*uy*uz - 1.*uz*uz + 
           0.3333333333333333*(ux*ux + 2.*uy + uy*uy + uz*(2. + uz))) - 0.03125*(-(nx*ux) - ny*(1 + uy) - nz*(1 + uz))*
          (2*chem*(uy + uz)*(uy + uz) - 0.3333333333333333*(-((rhoA - rhoB)*(uy + uz)*(uy + uz)) + 
             2*chem*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))); 
                feq17 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) - 
         0.0625*(Fz + Fx*ux + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy + 2.*uy*uz - 1.*uz*uz + 
           0.3333333333333333*(ux*ux - 2.*uy + uy*uy + uz*(2. + uz))) - 0.03125*(ny - nx*ux - ny*uy - nz*(1 + uz))*
          (2*chem*(uy - uz)*(uy - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(uy - uz)*(uy - uz)) + 
             2*chem*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))); 
                feq18 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) - 
         0.0625*(Fx*ux + Fy*(1 + uy) + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*uy*uy + 2.*uy*uz - 1.*uz*uz + 
           0.3333333333333333*(ux*ux + 2.*uy + uy*uy + (-2. + uz)*uz)) - 0.03125*(nz - nx*ux - ny*(1 + uy) - nz*uz)*
          (2*chem*(uy - uz)*(uy - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(uy - uz)*(uy - uz)) + 
             2*chem*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz))); 

            //------------------------------------------------- BCK collison ------------------------------------------------------------//
    		// q=0
    		dist[n] = m0 - (m0-feq0)/tau + 0.25*(2*(Fx*ux + Fy*uy + Fz*uz)*(-0.6666666666666666 + ux*ux + uy*uy + uz*uz) + 
         (mgx*ux + mgy*uy + mgz*uz)*(2*chem*(ux*ux + uy*uy + uz*uz) + 
            0.3333333333333333*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + uz*uz)))); 

    		// q = 1
    		dist[nr2] = m1 - (m1-feq1)/tau + 0.125*(2*(Fx*(-1 + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - ux*ux + 
            0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
         (mgx*(-1 + ux) + mgy*uy + mgz*uz)*(-2*chem*(ux*ux) + 
            0.3333333333333333*((-rhoA + rhoB)*(ux*ux) + 2*chem*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*uz))));

    		// q=2
    		dist[nr1] = m2 - (m2-feq2)/tau + 0.125*(2*(Fx + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - ux*ux + 
            0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*uz)) + 
         (mgx + mgx*ux + mgy*uy + mgz*uz)*(-2*chem*(ux*ux) + 
            0.3333333333333333*((-rhoA + rhoB)*(ux*ux) + 2*chem*(2*ux + ux*ux + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*uz))));

    		// q = 3
    		dist[nr4] = m3 - (m3-feq3)/tau + 0.125*(2*(Fx*ux + Fy*(-1 + uy) + Fz*uz)*(-0.2222222222222222 - uy*uy + 
            0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
         (mgx*ux + mgy*(-1 + uy) + mgz*uz)*(-2*chem*(uy*uy) + 
            0.3333333333333333*((-rhoA + rhoB)*(uy*uy) + 2*chem*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*uz))));

    		// q = 4
    		dist[nr3] = m4 - (m4-feq4)/tau + 0.125*(2*(Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - uy*uy + 
            0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*uz)) + 
         (mgy + mgx*ux + mgy*uy + mgz*uz)*(-2*chem*(uy*uy) + 
            0.3333333333333333*((-rhoA + rhoB)*(uy*uy) + 2*chem*(ux*ux + 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*uz))));

    		// q = 5
    		dist[nr6] = m5 - (m5-feq5)/tau + 0.125*(2*(Fx*ux + Fy*uy + Fz*(-1 + uz))*(-0.2222222222222222 - uz*uz + 
            0.3333333333333333*(ux*ux + uy*uy + (-2 + uz)*uz)) + 
         (mgx*ux + mgy*uy + mgz*(-1 + uz))*(-2*chem*(uz*uz) + 
            0.3333333333333333*((-rhoA + rhoB)*(uz*uz) + 2*chem*(ux*ux + uy*uy + (-2 + uz)*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + (-2 + uz)*uz))));

    		// q = 6
    		dist[nr5] = m6 - (m6-feq6)/tau + 0.125*(2*(Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - uz*uz + 
            0.3333333333333333*(ux*ux + uy*uy + uz*(2 + uz))) + 
         (mgz + mgx*ux + mgy*uy + mgz*uz)*(-2*chem*(uz*uz) + 
            0.3333333333333333*((-rhoA + rhoB)*(uz*uz) + 2*chem*(ux*ux + uy*uy + uz*(2 + uz))) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + uz*(2 + uz)))));

    		// q = 7
    		dist[nr8] = m7 - (m7-feq7)/tau + 0.0625*(-2*(Fx*(-1 + ux) + Fy*(-1 + uy) + Fz*uz)*
          (0.2222222222222222 + (ux + uy)*(ux + uy) - 
            0.3333333333333333*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
         (mgx*(-1 + ux) + mgy*(-1 + uy) + mgz*uz)*
          (-2*chem*((ux + uy)*(ux + uy)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux + uy)*(ux + uy))) + 2*chem*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz))));

    		// q = 8
    		dist[nr7] = m8 - (m8-feq8)/tau + 0.0625*(2*(Fx + Fy + Fx*ux + Fy*uy + Fz*uz)*
          (-0.2222222222222222 - (ux + uy)*(ux + uy) + 
            0.3333333333333333*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
         (mgx + mgy + mgx*ux + mgy*uy + mgz*uz)*
          (-2*chem*((ux + uy)*(ux + uy)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux + uy)*(ux + uy))) + 2*chem*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz))));

    		// q = 9
    		dist[nr10] = m9 - (m9-feq9)/tau + 0.0625*(2*(Fy + Fx*(-1 + ux) + Fy*uy + Fz*uz)*
          (-0.2222222222222222 - (ux - uy)*(ux - uy) + 
            0.3333333333333333*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
         (mgy + mgx*(-1 + ux) + mgy*uy + mgz*uz)*
          (-2*chem*((ux - uy)*(ux - uy)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux - uy)*(ux - uy))) + 2*chem*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz))));

    		// q = 10
    		dist[nr9] = m10 - (m10-feq10)/tau + 0.0625*(2*(Fx*(1 + ux) + Fy*(-1 + uy) + Fz*uz)*
          (-0.2222222222222222 - (ux - uy)*(ux - uy) + 
            0.3333333333333333*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
         (mgx*(1 + ux) + mgy*(-1 + uy) + mgz*uz)*
          (-2*chem*((ux - uy)*(ux - uy)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux - uy)*(ux - uy))) + 2*chem*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz))));

    		// q = 11
    		dist[nr12] = m11 - (m11-feq11)/tau + 0.0625*(-2*(Fx*(-1 + ux) + Fy*uy + Fz*(-1 + uz))*
          (0.2222222222222222 + (ux + uz)*(ux + uz) - 
            0.3333333333333333*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
         (mgx*(-1 + ux) + mgy*uy + mgz*(-1 + uz))*
          (-2*chem*((ux + uz)*(ux + uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux + uz)*(ux + uz))) + 2*chem*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz))));

    		// q = 12
    		dist[nr11] = m12 - (m12-feq12)/tau + 0.0625*(2*(Fx + Fz + Fx*ux + Fy*uy + Fz*uz)*
          (-0.2222222222222222 - (ux + uz)*(ux + uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))
           + (mgx + mgz + mgx*ux + mgy*uy + mgz*uz)*
          (-2*chem*((ux + uz)*(ux + uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux + uz)*(ux + uz))) + 2*chem*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))));

    		// q = 13
    		dist[nr14] = m13 - (m13-feq13)/tau + 0.0625*(2*(Fz + Fx*(-1 + ux) + Fy*uy + Fz*uz)*
          (-0.2222222222222222 - (ux - uz)*(ux - uz) + 
            0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
         (mgz + mgx*(-1 + ux) + mgy*uy + mgz*uz)*
          (-2*chem*((ux - uz)*(ux - uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux - uz)*(ux - uz))) + 2*chem*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*(2 + uz)))));

    		// q= 14
    		dist[nr13] = m14 - (m14-feq14)/tau + 0.0625*(2*(Fx*(1 + ux) + Fy*uy + Fz*(-1 + uz))*
          (-0.2222222222222222 - (ux - uz)*(ux - uz) + 
            0.3333333333333333*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
         (mgx*(1 + ux) + mgy*uy + mgz*(-1 + uz))*
          (-2*chem*((ux - uz)*(ux - uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux - uz)*(ux - uz))) + 2*chem*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz))));

    		// q = 15
    		dist[nr16] = m15 - (m15-feq15)/tau + 0.0625*(-2*(Fx*ux + Fy*(-1 + uy) + Fz*(-1 + uz))*
          (0.2222222222222222 + (uy + uz)*(uy + uz) - 0.3333333333333333*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))
           + (mgx*ux + mgy*(-1 + uy) + mgz*(-1 + uz))*
          (-2*chem*((uy + uz)*(uy + uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((uy + uz)*(uy + uz))) + 2*chem*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))));

    		// q = 16
    		dist[nr15] = m16 - (m16-feq16)/tau + 0.0625*(2*(Fy + Fz + Fx*ux + Fy*uy + Fz*uz)*
          (-0.2222222222222222 - (uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))
           + (mgy + mgz + mgx*ux + mgy*uy + mgz*uz)*
          (-2*chem*((uy + uz)*(uy + uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((uy + uz)*(uy + uz))) + 2*chem*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))));

    		// q = 17
    		dist[nr18] = m17 - (m17-feq17)/tau + 0.0625*(2*(Fz + Fx*ux + Fy*(-1 + uy) + Fz*uz)*
          (-0.2222222222222222 - (uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))
           + (mgz + mgx*ux + mgy*(-1 + uy) + mgz*uz)*
          (-2*chem*((uy - uz)*(uy - uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((uy - uz)*(uy - uz))) + 2*chem*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))));

    		// q = 18
    		dist[nr17] = m18 - (m18-feq18)/tau + 0.0625*(2*(Fx*ux + Fy*(1 + uy) + Fz*(-1 + uz))*
          (-0.2222222222222222 - (uy - uz)*(uy - uz) + 
            0.3333333333333333*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 
         (mgx*ux + mgy*(1 + uy) + mgz*(-1 + uz))*
          (-2*chem*((uy - uz)*(uy - uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((uy - uz)*(uy - uz))) + 2*chem*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz))));
            //----------------------------------------------------------------------------------------------------------------------------------------//

            //Update velocity on device
    		Vel[0*Np+n] = ux;
    		Vel[1*Np+n] = uy;
    		Vel[2*Np+n] = uz;
            //Update pressure on device
            Pressure[n] = p;
            //Update chemical potential on device
            mu_phi[n] = chem;
            //Update color gradient on device
    		ColorGrad[0*Np+n] = nx;
    		ColorGrad[1*Np+n] = ny;
    		ColorGrad[2*Np+n] = nz;
    	}
    }
}

__global__ void dvc_ScaLBL_D3Q19_AAeven_FreeLeeModel(int *Map, double *dist, double *Den,	double *Phi, double *mu_phi, double *Vel, double *Pressure, double *ColorGrad,
        double rhoA, double rhoB, double tauA, double tauB, double kappa, double beta, double W, double Fx, double Fy, double Fz, 
        int strideY, int strideZ, int start, int finish, int Np){

	int n,nn,nn2x,ijk;
    double ux,uy,uz;//fluid velocity 
    double p;//pressure
    double chem;//chemical potential
    double phi; //phase field
    double rho0;//fluid density
	// distribution functions
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
	double m0,m3,m5,m7;
	double mm1,mm2,mm4,mm6,mm8,mm9,mm10,mm11,mm12,mm13,mm14,mm15,mm16,mm17,mm18;
	double mm3,mm5,mm7;
    double feq0,feq1,feq2,feq3,feq4,feq5,feq6,feq7,feq8,feq9,feq10,feq11,feq12,feq13,feq14,feq15,feq16,feq17,feq18;
    double nx,ny,nz;//normal color gradient
    double mgx,mgy,mgz;//mixed gradient reaching secondary neighbor

	//double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18;
    //double h0,h1,h2,h3,h4,h5,h6;//distributions for LB phase field
	double tau;//position dependent LB relaxation time for fluid
    //double C,theta;
    //double M = 2.0/9.0*(tauM-0.5);//diffusivity (or mobility) for the phase field D3Q7

    //	for (int n=start; n<finish; n++){
    int S = Np/NBLOCKS/NTHREADS + 1;
    for (int s=0; s<S; s++){
    	//........Get 1-D index for this thread....................
    	n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

    	if ( n<finish ){
    		rho0 = Den[n];//load density

    		// Get the 1D index based on regular data layout
    		ijk = Map[n];
            phi = Phi[ijk];// load phase field

    		// local relaxation time
    		tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);

    		//					COMPUTE THE COLOR GRADIENT
    		//........................................................................
    		//.................Read Phase Indicator Values............................
    		//........................................................................
    		nn = ijk-1;							// neighbor index (get convention)
    		m1 = Phi[nn];						// get neighbor for phi - 1
    		//........................................................................
    		nn = ijk+1;							// neighbor index (get convention)
    		m2 = Phi[nn];						// get neighbor for phi - 2
    		//........................................................................
    		nn = ijk-strideY;							// neighbor index (get convention)
    		m3 = Phi[nn];					// get neighbor for phi - 3
    		//........................................................................
    		nn = ijk+strideY;							// neighbor index (get convention)
    		m4 = Phi[nn];					// get neighbor for phi - 4
    		//........................................................................
    		nn = ijk-strideZ;						// neighbor index (get convention)
    		m5 = Phi[nn];					// get neighbor for phi - 5
    		//........................................................................
    		nn = ijk+strideZ;						// neighbor index (get convention)
    		m6 = Phi[nn];					// get neighbor for phi - 6
    		//........................................................................
    		nn = ijk-strideY-1;						// neighbor index (get convention)
    		m7 = Phi[nn];					// get neighbor for phi - 7
    		//........................................................................
    		nn = ijk+strideY+1;						// neighbor index (get convention)
    		m8 = Phi[nn];					// get neighbor for phi - 8
    		//........................................................................
    		nn = ijk+strideY-1;						// neighbor index (get convention)
    		m9 = Phi[nn];					// get neighbor for phi - 9
    		//........................................................................
    		nn = ijk-strideY+1;						// neighbor index (get convention)
    		m10 = Phi[nn];					// get neighbor for phi - 10
    		//........................................................................
    		nn = ijk-strideZ-1;						// neighbor index (get convention)
    		m11 = Phi[nn];					// get neighbor for phi - 11
    		//........................................................................
    		nn = ijk+strideZ+1;						// neighbor index (get convention)
    		m12 = Phi[nn];					// get neighbor for phi - 12
    		//........................................................................
    		nn = ijk+strideZ-1;						// neighbor index (get convention)
    		m13 = Phi[nn];					// get neighbor for phi - 13
    		//........................................................................
    		nn = ijk-strideZ+1;						// neighbor index (get convention)
    		m14 = Phi[nn];					// get neighbor for phi - 14
    		//........................................................................
    		nn = ijk-strideZ-strideY;					// neighbor index (get convention)
    		m15 = Phi[nn];					// get neighbor for phi - 15
    		//........................................................................
    		nn = ijk+strideZ+strideY;					// neighbor index (get convention)
    		m16 = Phi[nn];					// get neighbor for phi - 16
    		//........................................................................
    		nn = ijk+strideZ-strideY;					// neighbor index (get convention)
    		m17 = Phi[nn];					// get neighbor for phi - 17
    		//........................................................................
    		nn = ijk-strideZ+strideY;					// neighbor index (get convention)
    		m18 = Phi[nn];					// get neighbor for phi - 18

            // compute mixed difference (Eq.30, A.Fukhari et al. JCP 315(2016) 434-457)
    		//........................................................................
    		nn2x = ijk-2;							// neighbor index (get convention)
    		mm1 = Phi[nn2x];						// get neighbor for phi - 1
            mm1 = 0.25*(-mm1+5.0*m1-3.0*phi-m2);
    		//........................................................................
    		nn2x = ijk+2;							// neighbor index (get convention)
    		mm2 = Phi[nn2x];						// get neighbor for phi - 2
            mm2 = 0.25*(-mm2+5.0*m2-3.0*phi-m1);
    		//........................................................................
    		nn2x = ijk-strideY*2;							// neighbor index (get convention)
    		mm3 = Phi[nn2x];					// get neighbor for phi - 3
            mm3 = 0.25*(-mm3+5.0*m3-3.0*phi-m4);
    		//........................................................................
    		nn2x = ijk+strideY*2;							// neighbor index (get convention)
    		mm4 = Phi[nn2x];					// get neighbor for phi - 4
            mm4 = 0.25*(-mm4+5.0*m4-3.0*phi-m3);
    		//........................................................................
    		nn2x = ijk-strideZ*2;						// neighbor index (get convention)
    		mm5 = Phi[nn2x];					// get neighbor for phi - 5
            mm5 = 0.25*(-mm5+5.0*m5-3.0*phi-m6);
    		//........................................................................
    		nn2x = ijk+strideZ*2;						// neighbor index (get convention)
    		mm6 = Phi[nn2x];					// get neighbor for phi - 6
            mm6 = 0.25*(-mm6+5.0*m6-3.0*phi-m5);
    		//........................................................................
    		nn2x = ijk-strideY*2-2;						// neighbor index (get convention)
    		mm7 = Phi[nn2x];					// get neighbor for phi - 7
            mm7 = 0.25*(-mm7+5.0*m7-3.0*phi-m8);
    		//........................................................................
    		nn2x = ijk+strideY*2+2;						// neighbor index (get convention)
    		mm8 = Phi[nn2x];					// get neighbor for phi - 8
            mm8 = 0.25*(-mm8+5.0*m8-3.0*phi-m7);
    		//........................................................................
    		nn2x = ijk+strideY*2-2;						// neighbor index (get convention)
    		mm9 = Phi[nn2x];					// get neighbor for phi - 9
            mm9 = 0.25*(-mm9+5.0*m9-3.0*phi-m10);
    		//........................................................................
    		nn2x = ijk-strideY*2+2;						// neighbor index (get convention)
    		mm10 = Phi[nn2x];					// get neighbor for phi - 10
            mm10 = 0.25*(-mm10+5.0*m10-3.0*phi-m9);
    		//........................................................................
    		nn2x = ijk-strideZ*2-2;						// neighbor index (get convention)
    		mm11 = Phi[nn2x];					// get neighbor for phi - 11
            mm11 = 0.25*(-mm11+5.0*m11-3.0*phi-m12);
    		//........................................................................
    		nn2x = ijk+strideZ*2+2;						// neighbor index (get convention)
    		mm12 = Phi[nn2x];					// get neighbor for phi - 12
            mm12 = 0.25*(-mm12+5.0*m12-3.0*phi-m11);
    		//........................................................................
    		nn2x = ijk+strideZ*2-2;						// neighbor index (get convention)
    		mm13 = Phi[nn2x];					// get neighbor for phi - 13
            mm13 = 0.25*(-mm13+5.0*m13-3.0*phi-m14);
    		//........................................................................
    		nn2x = ijk-strideZ*2+2;						// neighbor index (get convention)
    		mm14 = Phi[nn2x];					// get neighbor for phi - 14
            mm14 = 0.25*(-mm14+5.0*m14-3.0*phi-m13);
    		//........................................................................
    		nn2x = ijk-strideZ*2-strideY*2;					// neighbor index (get convention)
    		mm15 = Phi[nn2x];					// get neighbor for phi - 15
            mm15 = 0.25*(-mm15+5.0*m15-3.0*phi-m16);
    		//........................................................................
    		nn2x = ijk+strideZ*2+strideY*2;					// neighbor index (get convention)
    		mm16 = Phi[nn2x];					// get neighbor for phi - 16
            mm16 = 0.25*(-mm16+5.0*m16-3.0*phi-m15);
    		//........................................................................
    		nn2x = ijk+strideZ*2-strideY*2;					// neighbor index (get convention)
    		mm17 = Phi[nn2x];					// get neighbor for phi - 17
            mm17 = 0.25*(-mm17+5.0*m17-3.0*phi-m18);
    		//........................................................................
    		nn2x = ijk-strideZ*2+strideY*2;					// neighbor index (get convention)
    		mm18 = Phi[nn2x];					// get neighbor for phi - 18
            mm18 = 0.25*(-mm18+5.0*m18-3.0*phi-m17);


    		//............Compute the Color Gradient...................................
    		nx = -3.0*1.0/18.0*(m1-m2+0.5*(m7-m8+m9-m10+m11-m12+m13-m14));
    		ny = -3.0*1.0/18.0*(m3-m4+0.5*(m7-m8-m9+m10+m15-m16+m17-m18));
    		nz = -3.0*1.0/18.0*(m5-m6+0.5*(m11-m12-m13+m14+m15-m16-m17+m18));
    		//............Compute the Chemical Potential...............................
            chem = 2.0*3.0/18.0*(m1+m2+m3+m4+m5+m6-6*phi+0.5*(m7+m8+m9+m10+m11+m12+m13+m14+m15+m16+m17+m18-12*phi));//intermediate var, i.e. the laplacian
            chem = 4.0*beta*phi*(phi+1.0)*(phi-1.0)-kappa*chem;
    		//............Compute the Mixed Gradient...................................
    		mgx = -3.0*1.0/18.0*(mm1-mm2+0.5*(mm7-mm8+mm9-mm10+mm11-mm12+mm13-mm14));
    		mgy = -3.0*1.0/18.0*(mm3-mm4+0.5*(mm7-mm8-mm9+mm10+mm15-mm16+mm17-mm18));
    		mgz = -3.0*1.0/18.0*(mm5-mm6+0.5*(mm11-mm12-mm13+mm14+mm15-mm16-mm17+mm18));
    		
    		// q=0
    		m0 = dist[n];
    		// q=1
    		m1 = dist[2*Np+n]; 

            // q=2
    		m2 = dist[1*Np+n];  

    		// q=3
    		m3 = dist[4*Np+n];

    		// q = 4
    		m4 = dist[3*Np+n];

    		// q=5
    		m5 = dist[6*Np+n];

    		// q = 6
    		m6 = dist[5*Np+n];
    		
    		// q=7
    		m7 = dist[8*Np+n];

    		// q = 8
    		m8 = dist[7*Np+n];

    		// q=9
    		m9 = dist[10*Np+n];

    		// q = 10
    		m10 = dist[9*Np+n];

    		// q=11
    		m11 = dist[12*Np+n];

    		// q=12
    		m12 = dist[11*Np+n];

    		// q=13
    		m13 = dist[14*Np+n];

    		// q=14
    		m14 = dist[13*Np+n];

    		// q=15
    		m15 = dist[16*Np+n];

    		// q=16
    		m16 = dist[15*Np+n];

    		// q=17
    		m17 = dist[18*Np+n];

    		// q=18
    		m18 = dist[17*Np+n];

            //compute fluid velocity
            ux = 3.0/rho0*(m1-m2+m7-m8+m9-m10+m11-m12+m13-m14+0.5*(chem*nx+Fx)/3.0);
            uy = 3.0/rho0*(m3-m4+m7-m8-m9+m10+m15-m16+m17-m18+0.5*(chem*ny+Fy)/3.0);
            uz = 3.0/rho0*(m5-m6+m11-m12-m13+m14+m15-m16-m17+m18+0.5*(chem*nz+Fz)/3.0);

            //compute pressure
            p = (m0+m2+m1+m4+m3+m6+m5+m8+m7+m10+m9+m12+m11+m14+m13+m16+m15+m18+m17)
                      +0.5*(rhoA-rhoB)/2.0/3.0*(ux*nx+uy*ny+uz*nz);

            //compute equilibrium distributions
                feq0 = 0.3333333333333333*p - 0.25*(Fx*ux + Fy*uy + Fz*uz)*(-0.6666666666666666 + ux*ux + uy*uy + uz*uz) - 
         0.16666666666666666*rho0*(ux*ux + uy*uy + uz*uz) - 0.5*(-(nx*ux) - ny*uy - nz*uz)*
          (-0.08333333333333333*(rhoA - rhoB)*(ux*ux + uy*uy + uz*uz) + chem*(0.3333333333333333 - 0.5*(ux*ux + uy*uy + uz*uz)));
                feq1 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-ux*ux + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*uz)) - 
         0.125*(Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 
           0.3333333333333333*(-2.*ux + ux*ux + uy*uy + uz*uz)) - 0.0625*(nx - nx*ux - ny*uy - nz*uz)*
          (2*chem*ux*ux - 0.3333333333333333*((-rhoA + rhoB)*ux*ux + 2*chem*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
           0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*uz)));
                feq2 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-ux*ux + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*uz)) - 
         0.125*(Fx + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 
           0.3333333333333333*(2.*ux + ux*ux + uy*uy + uz*uz)) - 0.0625*(nx + nx*ux + ny*uy + nz*uz)*
          (-2.*chem*ux*ux + 0.1111111111111111*(-4.*chem + rhoB*(-2.*ux - 1.*ux*ux - 1.*uy*uy - 1.*uz*uz) + 
             rhoA*(2.*ux + ux*ux + uy*uy + uz*uz)) + 0.3333333333333333*((-1.*rhoA + rhoB)*ux*ux + 
             chem*(4.*ux + 2.*ux*ux + 2.*uy*uy + 2.*uz*uz)));
                feq3 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uy*uy + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*uz)) - 
         0.125*(Fx*ux + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy + 
           0.3333333333333333*(ux*ux - 2.*uy + uy*uy + uz*uz)) - 0.0625*(ny - nx*ux - ny*uy - nz*uz)*
          (2*chem*uy*uy - 0.3333333333333333*((-rhoA + rhoB)*uy*uy + 2*chem*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
           0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*uz)));
                feq4 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uy*uy + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*uz)) - 
         0.125*(Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy + 
           0.3333333333333333*(ux*ux + 2.*uy + uy*uy + uz*uz)) - 0.0625*(ny + nx*ux + ny*uy + nz*uz)*
          (-2.*chem*uy*uy + 0.1111111111111111*(-4.*chem + rhoB*(-1.*ux*ux - 2.*uy - 1.*uy*uy - 1.*uz*uz) + 
             rhoA*(ux*ux + 2.*uy + uy*uy + uz*uz)) + 0.3333333333333333*((-1.*rhoA + rhoB)*uy*uy + 
             chem*(2.*ux*ux + 4.*uy + 2.*uy*uy + 2.*uz*uz))); 
                feq5 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uz*uz + 0.3333333333333333*(ux*ux + uy*uy + (-2 + uz)*uz)) - 
         0.125*(Fx*ux + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*uz*uz + 
           0.3333333333333333*(ux*ux + uy*uy + (-2. + uz)*uz)) - 0.0625*(nx*ux + ny*uy + nz*(-1. + uz))*
          (-2.*chem*uz*uz + 0.1111111111111111*(-4.*chem + rhoB*(-1.*ux*ux - 1.*uy*uy + (2. - 1.*uz)*uz) + 
             rhoA*(ux*ux + uy*uy + (-2. + uz)*uz)) + 0.3333333333333333*((-1.*rhoA + rhoB)*uz*uz + 
             chem*(2.*ux*ux + 2.*uy*uy + uz*(-4. + 2.*uz)))); 
                feq6 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uz*uz + 0.3333333333333333*(ux*ux + uy*uy + uz*(2 + uz))) - 
         0.125*(Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*uz*uz + 
           0.3333333333333333*(ux*ux + uy*uy + uz*(2. + uz))) - 0.0625*(nz + nx*ux + ny*uy + nz*uz)*
          (-2.*chem*uz*uz + 0.1111111111111111*(-4.*chem + rhoB*(-1.*ux*ux - 1.*uy*uy + (-2. - 1.*uz)*uz) + 
             rhoA*(ux*ux + uy*uy + uz*(2. + uz))) + 0.3333333333333333*((-1.*rhoA + rhoB)*uz*uz + 
             chem*(2.*ux*ux + 2.*uy*uy + uz*(4. + 2.*uz)))); 
                feq7 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux + uy)*(ux + uy) + 0.3333333333333333*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) - 
         0.0625*(Fx*(-1. + ux) + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uy - 1.*uy*uy + 
           0.3333333333333333*(-2.*ux + ux*ux - 2.*uy + uy*uy + uz*uz)) - 0.03125*(nx + ny - nx*ux - ny*uy - nz*uz)*
          (2*chem*(ux + uy)*(ux + uy) + 0.3333333333333333*((rhoA - rhoB)*(ux + uy)*(ux + uy) - 2*chem*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
           0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz))); 
                feq8 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux + uy)*(ux + uy) + 0.3333333333333333*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) - 
         0.0625*(Fx + Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uy - 1.*uy*uy + 
           0.3333333333333333*(2.*ux + ux*ux + 2.*uy + uy*uy + uz*uz)) - 0.03125*(-(nx*(1 + ux)) - ny*(1 + uy) - nz*uz)*
          (2*chem*(ux + uy)*(ux + uy) - 0.3333333333333333*(-((rhoA - rhoB)*(ux + uy)*(ux + uy)) + 
             2*chem*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz))); 
                feq9 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux - uy)*(ux - uy) + 0.3333333333333333*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) - 
         0.0625*(Fy + Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uy - 1.*uy*uy + 
           0.3333333333333333*(-2.*ux + ux*ux + 2.*uy + uy*uy + uz*uz)) - 0.03125*(nx - nx*ux - ny*(1 + uy) - nz*uz)*
          (2*chem*(ux - uy)*(ux - uy) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uy)*(ux - uy)) + 
             2*chem*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz))); 
                feq10 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux - uy)*(ux - uy) + 0.3333333333333333*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) - 
         0.0625*(Fx*(1 + ux) + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uy - 1.*uy*uy + 
           0.3333333333333333*(2.*ux + ux*ux - 2.*uy + uy*uy + uz*uz)) - 0.03125*(ny - nx*(1 + ux) - ny*uy - nz*uz)*
          (2*chem*(ux - uy)*(ux - uy) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uy)*(ux - uy)) + 
             2*chem*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz))); 
                feq11 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux + uz)*(ux + uz) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) - 
         0.0625*(Fx*(-1. + ux) + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uz - 1.*uz*uz + 
           0.3333333333333333*(-2.*ux + ux*ux + uy*uy + (-2. + uz)*uz)) - 0.03125*(nx + nz - nx*ux - ny*uy - nz*uz)*
          (2*chem*(ux + uz)*(ux + uz) + 0.3333333333333333*((rhoA - rhoB)*(ux + uz)*(ux + uz) - 2*chem*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
           0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz))); 
                feq12 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux + uz)*(ux + uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) - 
         0.0625*(Fx + Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uz - 1.*uz*uz + 
           0.3333333333333333*(2.*ux + ux*ux + uy*uy + uz*(2. + uz))) - 0.03125*(-(nx*(1 + ux)) - ny*uy - nz*(1 + uz))*
          (2*chem*(ux + uz)*(ux + uz) - 0.3333333333333333*(-((rhoA - rhoB)*(ux + uz)*(ux + uz)) + 
             2*chem*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))); 
                feq13 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux - uz)*(ux - uz) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) - 
         0.0625*(Fz + Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uz - 1.*uz*uz + 
           0.3333333333333333*(-2.*ux + ux*ux + uy*uy + uz*(2. + uz))) - 0.03125*(nx - nx*ux - ny*uy - nz*(1 + uz))*
          (2*chem*(ux - uz)*(ux - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uz)*(ux - uz)) + 
             2*chem*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*(2 + uz)))); 
                feq14 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux - uz)*(ux - uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) - 
         0.0625*(Fx*(1 + ux) + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uz - 1.*uz*uz + 
           0.3333333333333333*(2.*ux + ux*ux + uy*uy + (-2. + uz)*uz)) - 0.03125*(nz - nx*(1 + ux) - ny*uy - nz*uz)*
          (2*chem*(ux - uz)*(ux - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uz)*(ux - uz)) + 
             2*chem*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz))); 
                feq15 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) - 
         0.0625*(Fx*ux + Fy*(-1. + uy) + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*uy*uy - 2.*uy*uz - 1.*uz*uz + 
           0.3333333333333333*(ux*ux - 2.*uy + uy*uy + (-2. + uz)*uz)) - 0.03125*(ny + nz - nx*ux - ny*uy - nz*uz)*
          (2*chem*(uy + uz)*(uy + uz) + 0.3333333333333333*((rhoA - rhoB)*(uy + uz)*(uy + uz) - 2*chem*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) + 
           0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))); 
                feq16 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) - 
         0.0625*(Fy + Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy - 2.*uy*uz - 1.*uz*uz + 
           0.3333333333333333*(ux*ux + 2.*uy + uy*uy + uz*(2. + uz))) - 0.03125*(-(nx*ux) - ny*(1 + uy) - nz*(1 + uz))*
          (2*chem*(uy + uz)*(uy + uz) - 0.3333333333333333*(-((rhoA - rhoB)*(uy + uz)*(uy + uz)) + 
             2*chem*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))); 
                feq17 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) - 
         0.0625*(Fz + Fx*ux + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy + 2.*uy*uz - 1.*uz*uz + 
           0.3333333333333333*(ux*ux - 2.*uy + uy*uy + uz*(2. + uz))) - 0.03125*(ny - nx*ux - ny*uy - nz*(1 + uz))*
          (2*chem*(uy - uz)*(uy - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(uy - uz)*(uy - uz)) + 
             2*chem*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))); 
                feq18 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) - 
         0.0625*(Fx*ux + Fy*(1 + uy) + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*uy*uy + 2.*uy*uz - 1.*uz*uz + 
           0.3333333333333333*(ux*ux + 2.*uy + uy*uy + (-2. + uz)*uz)) - 0.03125*(nz - nx*ux - ny*(1 + uy) - nz*uz)*
          (2*chem*(uy - uz)*(uy - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(uy - uz)*(uy - uz)) + 
             2*chem*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz))); 

            //------------------------------------------------- BCK collison ------------------------------------------------------------//
    		// q=0
    		dist[n] = m0 - (m0-feq0)/tau + 0.25*(2*(Fx*ux + Fy*uy + Fz*uz)*(-0.6666666666666666 + ux*ux + uy*uy + uz*uz) + 
         (mgx*ux + mgy*uy + mgz*uz)*(2*chem*(ux*ux + uy*uy + uz*uz) + 
            0.3333333333333333*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + uz*uz))));

    		// q = 1
    		dist[1*Np+n] = m1 - (m1-feq1)/tau + 0.125*(2*(Fx*(-1 + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - ux*ux + 
            0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
         (mgx*(-1 + ux) + mgy*uy + mgz*uz)*(-2*chem*(ux*ux) + 
            0.3333333333333333*((-rhoA + rhoB)*(ux*ux) + 2*chem*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*uz))));

    		// q=2
    		dist[2*Np+n] = m2 - (m2-feq2)/tau + 0.125*(2*(Fx + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - ux*ux + 
            0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*uz)) + 
         (mgx + mgx*ux + mgy*uy + mgz*uz)*(-2*chem*(ux*ux) + 
            0.3333333333333333*((-rhoA + rhoB)*(ux*ux) + 2*chem*(2*ux + ux*ux + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*uz))));

    		// q = 3
    		dist[3*Np+n] = m3 - (m3-feq3)/tau + 0.125*(2*(Fx*ux + Fy*(-1 + uy) + Fz*uz)*(-0.2222222222222222 - uy*uy + 
            0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
         (mgx*ux + mgy*(-1 + uy) + mgz*uz)*(-2*chem*(uy*uy) + 
            0.3333333333333333*((-rhoA + rhoB)*(uy*uy) + 2*chem*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*uz))));

    		// q = 4
    		dist[4*Np+n] = m4 - (m4-feq4)/tau + 0.125*(2*(Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - uy*uy + 
            0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*uz)) + 
         (mgy + mgx*ux + mgy*uy + mgz*uz)*(-2*chem*(uy*uy) + 
            0.3333333333333333*((-rhoA + rhoB)*(uy*uy) + 2*chem*(ux*ux + 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*uz))));

    		// q = 5
    		dist[5*Np+n] = m5 - (m5-feq5)/tau + 0.125*(2*(Fx*ux + Fy*uy + Fz*(-1 + uz))*(-0.2222222222222222 - uz*uz + 
            0.3333333333333333*(ux*ux + uy*uy + (-2 + uz)*uz)) + 
         (mgx*ux + mgy*uy + mgz*(-1 + uz))*(-2*chem*(uz*uz) + 
            0.3333333333333333*((-rhoA + rhoB)*(uz*uz) + 2*chem*(ux*ux + uy*uy + (-2 + uz)*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + (-2 + uz)*uz))));

    		// q = 6
    		dist[6*Np+n] = m6 - (m6-feq6)/tau + 0.125*(2*(Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - uz*uz + 
            0.3333333333333333*(ux*ux + uy*uy + uz*(2 + uz))) + 
         (mgz + mgx*ux + mgy*uy + mgz*uz)*(-2*chem*(uz*uz) + 
            0.3333333333333333*((-rhoA + rhoB)*(uz*uz) + 2*chem*(ux*ux + uy*uy + uz*(2 + uz))) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + uz*(2 + uz)))));

    		// q = 7
    		dist[7*Np+n] = m7 - (m7-feq7)/tau + 0.0625*(-2*(Fx*(-1 + ux) + Fy*(-1 + uy) + Fz*uz)*
          (0.2222222222222222 + (ux + uy)*(ux + uy) - 
            0.3333333333333333*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
         (mgx*(-1 + ux) + mgy*(-1 + uy) + mgz*uz)*
          (-2*chem*((ux + uy)*(ux + uy)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux + uy)*(ux + uy))) + 2*chem*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz))));

    		// q = 8
    		dist[8*Np+n] = m8 - (m8-feq8)/tau + 0.0625*(2*(Fx + Fy + Fx*ux + Fy*uy + Fz*uz)*
          (-0.2222222222222222 - (ux + uy)*(ux + uy) + 
            0.3333333333333333*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
         (mgx + mgy + mgx*ux + mgy*uy + mgz*uz)*
          (-2*chem*((ux + uy)*(ux + uy)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux + uy)*(ux + uy))) + 2*chem*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz))));

    		// q = 9
    		dist[9*Np+n] = m9 - (m9-feq9)/tau + 0.0625*(2*(Fy + Fx*(-1 + ux) + Fy*uy + Fz*uz)*
          (-0.2222222222222222 - (ux - uy)*(ux - uy) + 
            0.3333333333333333*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
         (mgy + mgx*(-1 + ux) + mgy*uy + mgz*uz)*
          (-2*chem*((ux - uy)*(ux - uy)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux - uy)*(ux - uy))) + 2*chem*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz))));

    		// q = 10
    		dist[10*Np+n] = m10 - (m10-feq10)/tau + 0.0625*(2*(Fx*(1 + ux) + Fy*(-1 + uy) + Fz*uz)*
          (-0.2222222222222222 - (ux - uy)*(ux - uy) + 
            0.3333333333333333*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
         (mgx*(1 + ux) + mgy*(-1 + uy) + mgz*uz)*
          (-2*chem*((ux - uy)*(ux - uy)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux - uy)*(ux - uy))) + 2*chem*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz))));

    		// q = 11
    		dist[11*Np+n] = m11 - (m11-feq11)/tau + 0.0625*(-2*(Fx*(-1 + ux) + Fy*uy + Fz*(-1 + uz))*
          (0.2222222222222222 + (ux + uz)*(ux + uz) - 
            0.3333333333333333*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
         (mgx*(-1 + ux) + mgy*uy + mgz*(-1 + uz))*
          (-2*chem*((ux + uz)*(ux + uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux + uz)*(ux + uz))) + 2*chem*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz))));

    		// q = 12
    		dist[12*Np+n] = m12 - (m12-feq12)/tau + 0.0625*(2*(Fx + Fz + Fx*ux + Fy*uy + Fz*uz)*
          (-0.2222222222222222 - (ux + uz)*(ux + uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))
           + (mgx + mgz + mgx*ux + mgy*uy + mgz*uz)*
          (-2*chem*((ux + uz)*(ux + uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux + uz)*(ux + uz))) + 2*chem*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))));

    		// q = 13
    		dist[13*Np+n] = m13 - (m13-feq13)/tau + 0.0625*(2*(Fz + Fx*(-1 + ux) + Fy*uy + Fz*uz)*
          (-0.2222222222222222 - (ux - uz)*(ux - uz) + 
            0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
         (mgz + mgx*(-1 + ux) + mgy*uy + mgz*uz)*
          (-2*chem*((ux - uz)*(ux - uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux - uz)*(ux - uz))) + 2*chem*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*(2 + uz)))));

    		// q= 14
    		dist[14*Np+n] = m14 - (m14-feq14)/tau + 0.0625*(2*(Fx*(1 + ux) + Fy*uy + Fz*(-1 + uz))*
          (-0.2222222222222222 - (ux - uz)*(ux - uz) + 
            0.3333333333333333*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
         (mgx*(1 + ux) + mgy*uy + mgz*(-1 + uz))*
          (-2*chem*((ux - uz)*(ux - uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux - uz)*(ux - uz))) + 2*chem*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz))));

    		// q = 15
    		dist[15*Np+n] = m15 - (m15-feq15)/tau + 0.0625*(-2*(Fx*ux + Fy*(-1 + uy) + Fz*(-1 + uz))*
          (0.2222222222222222 + (uy + uz)*(uy + uz) - 0.3333333333333333*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))
           + (mgx*ux + mgy*(-1 + uy) + mgz*(-1 + uz))*
          (-2*chem*((uy + uz)*(uy + uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((uy + uz)*(uy + uz))) + 2*chem*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))));

    		// q = 16
    		dist[16*Np+n] = m16 - (m16-feq16)/tau + 0.0625*(2*(Fy + Fz + Fx*ux + Fy*uy + Fz*uz)*
          (-0.2222222222222222 - (uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))
           + (mgy + mgz + mgx*ux + mgy*uy + mgz*uz)*
          (-2*chem*((uy + uz)*(uy + uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((uy + uz)*(uy + uz))) + 2*chem*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))));

    		// q = 17
    		dist[17*Np+n] = m17 - (m17-feq17)/tau + 0.0625*(2*(Fz + Fx*ux + Fy*(-1 + uy) + Fz*uz)*
          (-0.2222222222222222 - (uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))
           + (mgz + mgx*ux + mgy*(-1 + uy) + mgz*uz)*
          (-2*chem*((uy - uz)*(uy - uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((uy - uz)*(uy - uz))) + 2*chem*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))));

    		// q = 18
    		dist[18*Np+n] = m18 - (m18-feq18)/tau + 0.0625*(2*(Fx*ux + Fy*(1 + uy) + Fz*(-1 + uz))*
          (-0.2222222222222222 - (uy - uz)*(uy - uz) + 
            0.3333333333333333*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 
         (mgx*ux + mgy*(1 + uy) + mgz*(-1 + uz))*
          (-2*chem*((uy - uz)*(uy - uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((uy - uz)*(uy - uz))) + 2*chem*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz))));
            //----------------------------------------------------------------------------------------------------------------------------------------//

            //Update velocity on device
    		Vel[0*Np+n] = ux;
    		Vel[1*Np+n] = uy;
    		Vel[2*Np+n] = uz;
            //Update pressure on device
            Pressure[n] = p;
            //Update chemical potential on device
            mu_phi[n] = chem;
            //Update color gradient on device
    		ColorGrad[0*Np+n] = nx;
    		ColorGrad[1*Np+n] = ny;
    		ColorGrad[2*Np+n] = nz;

    	}
	}
}

__global__ void dvc_ScaLBL_D3Q19_AAodd_FreeLeeModel_Combined(int *neighborList, int *Map, double *dist, double *hq, double *Den,	double *Phi, double *mu_phi, double *Vel, double *Pressure, double *ColorGrad, 
        double rhoA, double rhoB, double tauA, double tauB, double tauM, double kappa, double beta, double W, double Fx, double Fy, double Fz, 
        int strideY, int strideZ, int start, int finish, int Np){

	int n,nn,nn2x,ijk;
	int nr1,nr2,nr3,nr4,nr5,nr6,nr7,nr8,nr9,nr10,nr11,nr12,nr13,nr14,nr15,nr16,nr17,nr18;
    double ux,uy,uz;//fluid velocity 
    double p;//pressure
    double chem;//chemical potential
    double phi; //phase field
    double rho0;//fluid density
	// distribution functions
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
	double m0,m3,m5,m7;
	//double mm1,mm2,mm4,mm6,mm8,mm9,mm10,mm11,mm12,mm13,mm14,mm15,mm16,mm17,mm18;
	//double mm3,mm5,mm7;
    double feq0,feq1,feq2,feq3,feq4,feq5,feq6,feq7,feq8,feq9,feq10,feq11,feq12,feq13,feq14,feq15,feq16,feq17,feq18;
    double nx,ny,nz;//normal color gradient
    double mgx,mgy,mgz;//mixed gradient reaching secondary neighbor
    double dirGradC1,dirGradC2,dirGradC3,dirGradC4,dirGradC5,dirGradC6,dirGradC7,dirGradC8,dirGradC9,dirGradC10,dirGradC11,dirGradC12;
    double dirGradC13,dirGradC14,dirGradC15,dirGradC16,dirGradC17,dirGradC18;
    double dirGradM1,dirGradM2,dirGradM3,dirGradM4,dirGradM5,dirGradM6,dirGradM7,dirGradM8,dirGradM9,dirGradM10,dirGradM11,dirGradM12;
    double dirGradM13,dirGradM14,dirGradM15,dirGradM16,dirGradM17,dirGradM18;

    double h0,h1,h2,h3,h4,h5,h6;//distributions for LB phase field
	double tau;//position dependent LB relaxation time for fluid
    double C,theta;
     double M = 2.0/9.0*(tauM-0.5);//diffusivity (or mobility) for the phase field D3Q7
    double phi_temp;

    int S = Np/NBLOCKS/NTHREADS + 1;
    for (int s=0; s<S; s++){
    	//........Get 1-D index for this thread....................
    	n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

    	if ( n<finish ){
    		rho0 = Den[n];//load density

    		// Get the 1D index based on regular data layout
    		ijk = Map[n];
            phi = Phi[ijk];// load phase field
            phi_temp = phi;
            if (phi>1.f) phi_temp=1.0;
            if (phi<-1.f) phi_temp=-1.0;

    		// local relaxation time
    		tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);

    		//					COMPUTE THE COLOR GRADIENT
    		//........................................................................
    		//.................Read Phase Indicator Values............................
    		//........................................................................
    		nn = ijk-1;							// neighbor index (get convention)
    		m2 = Phi[nn];						// get neighbor for phi - 1
    		//........................................................................
    		nn = ijk+1;							// neighbor index (get convention)
    		m1 = Phi[nn];						// get neighbor for phi - 2
            dirGradC1 = 0.5*(m1-m2);
            dirGradC2 = 0.5*(m2-m1);
    		//........................................................................
    		nn = ijk-strideY;							// neighbor index (get convention)
    		m4 = Phi[nn];					// get neighbor for phi - 3
    		//........................................................................
    		nn = ijk+strideY;							// neighbor index (get convention)
    		m3 = Phi[nn];					// get neighbor for phi - 4
            dirGradC3 = 0.5*(m3-m4);
            dirGradC4 = 0.5*(m4-m3);
    		//........................................................................
    		nn = ijk-strideZ;						// neighbor index (get convention)
    		m6 = Phi[nn];					// get neighbor for phi - 5
    		//........................................................................
    		nn = ijk+strideZ;						// neighbor index (get convention)
    		m5 = Phi[nn];					// get neighbor for phi - 6
            dirGradC5 = 0.5*(m5-m6);
            dirGradC6 = 0.5*(m6-m5);
    		//........................................................................
    		nn = ijk-strideY-1;						// neighbor index (get convention)
    		m8 = Phi[nn];					// get neighbor for phi - 7
    		//........................................................................
    		nn = ijk+strideY+1;						// neighbor index (get convention)
    		m7 = Phi[nn];					// get neighbor for phi - 8
            dirGradC7 = 0.5*(m7-m8);
            dirGradC8 = 0.5*(m8-m7);
    		//........................................................................
    		nn = ijk+strideY-1;						// neighbor index (get convention)
    		m10 = Phi[nn];					// get neighbor for phi - 9
    		//........................................................................
    		nn = ijk-strideY+1;						// neighbor index (get convention)
    		m9 = Phi[nn];					// get neighbor for phi - 10
            dirGradC9  = 0.5*(m9-m10);
            dirGradC10 = 0.5*(m10-m9);
    		//........................................................................
    		nn = ijk-strideZ-1;						// neighbor index (get convention)
    		m12 = Phi[nn];					// get neighbor for phi - 11
    		//........................................................................
    		nn = ijk+strideZ+1;						// neighbor index (get convention)
    		m11 = Phi[nn];					// get neighbor for phi - 12
            dirGradC11 = 0.5*(m11-m12);
            dirGradC12 = 0.5*(m12-m11);
    		//........................................................................
    		nn = ijk+strideZ-1;						// neighbor index (get convention)
    		m14 = Phi[nn];					// get neighbor for phi - 13
    		//........................................................................
    		nn = ijk-strideZ+1;						// neighbor index (get convention)
    		m13 = Phi[nn];					// get neighbor for phi - 14
            dirGradC13 = 0.5*(m13-m14);
            dirGradC14 = 0.5*(m14-m13);
    		//........................................................................
    		nn = ijk-strideZ-strideY;					// neighbor index (get convention)
    		m16 = Phi[nn];					// get neighbor for phi - 15
    		//........................................................................
    		nn = ijk+strideZ+strideY;					// neighbor index (get convention)
    		m15 = Phi[nn];					// get neighbor for phi - 16
            dirGradC15 = 0.5*(m15-m16);
            dirGradC16 = 0.5*(m16-m15);
    		//........................................................................
    		nn = ijk+strideZ-strideY;					// neighbor index (get convention)
    		m18 = Phi[nn];					// get neighbor for phi - 17
    		//........................................................................
    		nn = ijk-strideZ+strideY;					// neighbor index (get convention)
    		m17 = Phi[nn];					// get neighbor for phi - 18
            dirGradC17 = 0.5*(m17-m18);
            dirGradC18 = 0.5*(m18-m17);

            // compute mixed difference (Eq.30, A.Fukhari et al. JCP 315(2016) 434-457)
    		//........................................................................
    		nn2x = ijk+2;							// neighbor index (get convention)
    		dirGradM1 = Phi[nn2x];						// get neighbor for phi - 1
            dirGradM1 = 0.25*(-dirGradM1+5.0*m1-3.0*phi-m2);
    		//........................................................................
    		nn2x = ijk-2;							// neighbor index (get convention)
    		dirGradM2 = Phi[nn2x];						// get neighbor for phi - 2
            dirGradM2 = 0.25*(-dirGradM2+5.0*m2-3.0*phi-m1);
    		//........................................................................
    		nn2x = ijk+strideY*2;							// neighbor index (get convention)
    		dirGradM3 = Phi[nn2x];					// get neighbor for phi - 3
            dirGradM3 = 0.25*(-dirGradM3+5.0*m3-3.0*phi-m4);
    		//........................................................................
    		nn2x = ijk-strideY*2;							// neighbor index (get convention)
    		dirGradM4 = Phi[nn2x];					// get neighbor for phi - 4
            dirGradM4 = 0.25*(-dirGradM4+5.0*m4-3.0*phi-m3);
    		//........................................................................
    		nn2x = ijk+strideZ*2;						// neighbor index (get convention)
    		dirGradM5 = Phi[nn2x];					// get neighbor for phi - 5
            dirGradM5 = 0.25*(-dirGradM5+5.0*m5-3.0*phi-m6);
    		//........................................................................
    		nn2x = ijk-strideZ*2;						// neighbor index (get convention)
    		dirGradM6 = Phi[nn2x];					// get neighbor for phi - 6
            dirGradM6 = 0.25*(-dirGradM6+5.0*m6-3.0*phi-m5);
    		//........................................................................
    		nn2x = ijk+strideY*2+2;						// neighbor index (get convention)
    		dirGradM7 = Phi[nn2x];					// get neighbor for phi - 7
            dirGradM7 = 0.25*(-dirGradM7+5.0*m7-3.0*phi-m8);
    		//........................................................................
    		nn2x = ijk-strideY*2-2;						// neighbor index (get convention)
    		dirGradM8 = Phi[nn2x];					// get neighbor for phi - 8
            dirGradM8 = 0.25*(-dirGradM8+5.0*m8-3.0*phi-m7);
    		//........................................................................
    		nn2x = ijk-strideY*2+2;						// neighbor index (get convention)
    		dirGradM9 = Phi[nn2x];					// get neighbor for phi - 9
            dirGradM9 = 0.25*(-dirGradM9+5.0*m9-3.0*phi-m10);
    		//........................................................................
    		nn2x = ijk+strideY*2-2;						// neighbor index (get convention)
    		dirGradM10 = Phi[nn2x];					// get neighbor for phi - 10
            dirGradM10 = 0.25*(-dirGradM10+5.0*m10-3.0*phi-m9);
    		//........................................................................
    		nn2x = ijk+strideZ*2+2;						// neighbor index (get convention)
    		dirGradM11 = Phi[nn2x];					// get neighbor for phi - 11
            dirGradM11 = 0.25*(-dirGradM11+5.0*m11-3.0*phi-m12);
    		//........................................................................
    		nn2x = ijk-strideZ*2-2;						// neighbor index (get convention)
    		dirGradM12 = Phi[nn2x];					// get neighbor for phi - 12
            dirGradM12 = 0.25*(-dirGradM12+5.0*m12-3.0*phi-m11);
    		//........................................................................
    		nn2x = ijk-strideZ*2+2;						// neighbor index (get convention)
    		dirGradM13 = Phi[nn2x];					// get neighbor for phi - 13
            dirGradM13 = 0.25*(-dirGradM13+5.0*m13-3.0*phi-m14);
    		//........................................................................
    		nn2x = ijk+strideZ*2-2;						// neighbor index (get convention)
    		dirGradM14 = Phi[nn2x];					// get neighbor for phi - 14
            dirGradM14 = 0.25*(-dirGradM14+5.0*m14-3.0*phi-m13);
    		//........................................................................
    		nn2x = ijk+strideZ*2+strideY*2;					// neighbor index (get convention)
    		dirGradM15 = Phi[nn2x];					// get neighbor for phi - 15
            dirGradM15 = 0.25*(-dirGradM15+5.0*m15-3.0*phi-m16);
    		//........................................................................
    		nn2x = ijk-strideZ*2-strideY*2;					// neighbor index (get convention)
    		dirGradM16 = Phi[nn2x];					// get neighbor for phi - 16
            dirGradM16 = 0.25*(-dirGradM16+5.0*m16-3.0*phi-m15);
    		//........................................................................
    		nn2x = ijk-strideZ*2+strideY*2;					// neighbor index (get convention)
    		dirGradM17 = Phi[nn2x];					// get neighbor for phi - 17
            dirGradM17 = 0.25*(-dirGradM17+5.0*m17-3.0*phi-m18);
    		//........................................................................
    		nn2x = ijk+strideZ*2-strideY*2;					// neighbor index (get convention)
    		dirGradM18 = Phi[nn2x];					// get neighbor for phi - 18
            dirGradM18 = 0.25*(-dirGradM18+5.0*m18-3.0*phi-m17);


    		//............Compute the Color Gradient...................................
    		nx = 3.0*1.0/18.0*(dirGradC1-dirGradC2+0.5*(dirGradC7-dirGradC8+dirGradC9-dirGradC10+dirGradC11-dirGradC12+dirGradC13-dirGradC14));
    		ny = 3.0*1.0/18.0*(dirGradC3-dirGradC4+0.5*(dirGradC7-dirGradC8-dirGradC9+dirGradC10+dirGradC15-dirGradC16+dirGradC17-dirGradC18));
    		nz = 3.0*1.0/18.0*(dirGradC5-dirGradC6+0.5*(dirGradC11-dirGradC12-dirGradC13+dirGradC14+dirGradC15-dirGradC16-dirGradC17+dirGradC18));
    		//............Compute the Chemical Potential...............................
            //chem = 2.0*3.0/18.0*(m1+m2+m3+m4+m5+m6-6*phi+0.5*(m7+m8+m9+m10+m11+m12+m13+m14+m15+m16+m17+m18-12*phi));//intermediate var, i.e. the laplacian
            //chem = 4.0*beta*phi*(phi+1.0)*(phi-1.0)-kappa*chem;
            chem = 2.0*3.0/18.0*(m1+m2+m3+m4+m5+m6-6*phi_temp+0.5*(m7+m8+m9+m10+m11+m12+m13+m14+m15+m16+m17+m18-12*phi_temp));//intermediate var, i.e. the laplacian
            chem = 4.0*beta*phi_temp*(phi_temp+1.0)*(phi_temp-1.0)-kappa*chem;
    		//............Compute the Mixed Gradient...................................
    		mgx = 3.0*1.0/18.0*(dirGradM1-dirGradM2+0.5*(dirGradM7-dirGradM8+dirGradM9-dirGradM10+dirGradM11-dirGradM12+dirGradM13-dirGradM14));
    		mgy = 3.0*1.0/18.0*(dirGradM3-dirGradM4+0.5*(dirGradM7-dirGradM8-dirGradM9+dirGradM10+dirGradM15-dirGradM16+dirGradM17-dirGradM18));
    		mgz = 3.0*1.0/18.0*(dirGradM5-dirGradM6+0.5*(dirGradM11-dirGradM12-dirGradM13+dirGradM14+dirGradM15-dirGradM16-dirGradM17+dirGradM18));
    		
            //de-noise color gradient and mixed gradient
            C = sqrt(nx*nx+ny*ny+nz*nz);
            if (C<1.0e-12) nx=ny=nz=0.0;
            double mg_mag = sqrt(mgx*mgx+mgy*mgy+mgz*mgz);
            if (mg_mag<1.0e-12) mgx=mgy=mgz=0.0;
            //maybe you can also de-noise chemical potential ? within the bulk phase chem should be zero
            if (fabs(chem)<1.0e-12) chem=0.0;

    		// q=0
    		m0 = dist[n];
    		// q=1
    		nr1 = neighborList[n]; // neighbor 2 ( > 10Np => odd part of dist)
    		m1 = dist[nr1]; // reading the f1 data into register fq

    		nr2 = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
    		m2 = dist[nr2];  // reading the f2 data into register fq

    		// q=3
    		nr3 = neighborList[n+2*Np]; // neighbor 4
    		m3 = dist[nr3];

    		// q = 4
    		nr4 = neighborList[n+3*Np]; // neighbor 3
    		m4 = dist[nr4];

    		// q=5
    		nr5 = neighborList[n+4*Np];
    		m5 = dist[nr5];

    		// q = 6
    		nr6 = neighborList[n+5*Np];
    		m6 = dist[nr6];
    		
    		// q=7
    		nr7 = neighborList[n+6*Np];
    		m7 = dist[nr7];

    		// q = 8
    		nr8 = neighborList[n+7*Np];
    		m8 = dist[nr8];

    		// q=9
    		nr9 = neighborList[n+8*Np];
    		m9 = dist[nr9];

    		// q = 10
    		nr10 = neighborList[n+9*Np];
    		m10 = dist[nr10];

    		// q=11
    		nr11 = neighborList[n+10*Np];
    		m11 = dist[nr11];

    		// q=12
    		nr12 = neighborList[n+11*Np];
    		m12 = dist[nr12];

    		// q=13
    		nr13 = neighborList[n+12*Np];
    		m13 = dist[nr13];

    		// q=14
    		nr14 = neighborList[n+13*Np];
    		m14 = dist[nr14];

    		// q=15
    		nr15 = neighborList[n+14*Np];
    		m15 = dist[nr15];

    		// q=16
    		nr16 = neighborList[n+15*Np];
    		m16 = dist[nr16];

    		// q=17
    		nr17 = neighborList[n+16*Np];
    		m17 = dist[nr17];

    		// q=18
    		nr18 = neighborList[n+17*Np];
    		m18 = dist[nr18];

            //compute fluid velocity
            ux = 3.0/rho0*(m1-m2+m7-m8+m9-m10+m11-m12+m13-m14+0.5*(chem*nx+Fx)/3.0);
            uy = 3.0/rho0*(m3-m4+m7-m8-m9+m10+m15-m16+m17-m18+0.5*(chem*ny+Fy)/3.0);
            uz = 3.0/rho0*(m5-m6+m11-m12-m13+m14+m15-m16-m17+m18+0.5*(chem*nz+Fz)/3.0);
            //compute pressure
            p = (m0+m2+m1+m4+m3+m6+m5+m8+m7+m10+m9+m12+m11+m14+m13+m16+m15+m18+m17)
                      +0.5*(rhoA-rhoB)/2.0/3.0*(ux*nx+uy*ny+uz*nz);

            //compute equilibrium distributions
                feq0 = 0.3333333333333333*p - 0.25*(Fx*ux + Fy*uy + Fz*uz)*(-0.6666666666666666 + ux*ux + uy*uy + uz*uz) - 0.16666666666666666*rho0*(ux*ux + uy*uy + uz*uz) - 
   0.5*(-(nx*ux) - ny*uy - nz*uz)*(-0.08333333333333333*(rhoA - rhoB)*(ux*ux + uy*uy + uz*uz) + chem*(0.3333333333333333 - 0.5*(ux*ux + uy*uy + uz*uz))); 
                feq1 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-(ux*ux) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*uz)) - 
   0.125*(Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(ux*ux) + 0.3333333333333333*(-2.*ux + ux*ux + uy*uy + uz*uz)) - 
   0.0625*(dirGradC1 - nx*ux - ny*uy - nz*uz)*(2*chem*(ux*ux) - 0.3333333333333333*((-rhoA + rhoB)*(ux*ux) + 2*chem*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*uz))); 
                feq2 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-(ux*ux) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*uz)) - 
   0.125*(Fx + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(ux*ux) + 0.3333333333333333*(2.*ux + ux*ux + uy*uy + uz*uz)) - 
   0.0625*(dirGradC2 - nx*ux - ny*uy - nz*uz)*(2*chem*(ux*ux) - 0.3333333333333333*((-rhoA + rhoB)*(ux*ux) + 2*chem*(2*ux + ux*ux + uy*uy + uz*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*uz))); 
                feq3 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-(uy*uy) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*uz)) - 
   0.125*(Fx*ux + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*(uy*uy) + 0.3333333333333333*(ux*ux - 2.*uy + uy*uy + uz*uz)) - 
   0.0625*(dirGradC3 - nx*ux - ny*uy - nz*uz)*(2*chem*(uy*uy) - 0.3333333333333333*((-rhoA + rhoB)*(uy*uy) + 2*chem*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*uz)));  
                feq4 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-(uy*uy) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*uz)) - 
   0.125*(Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(uy*uy) + 0.3333333333333333*(ux*ux + 2.*uy + uy*uy + uz*uz)) - 
   0.0625*(dirGradC4 - nx*ux - ny*uy - nz*uz)*(2*chem*(uy*uy) - 0.3333333333333333*((-rhoA + rhoB)*(uy*uy) + 2*chem*(ux*ux + 2*uy + uy*uy + uz*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*uz)));  
                feq5 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-(uz*uz) + 0.3333333333333333*(ux*ux + uy*uy + (-2 + uz)*uz)) - 
   0.125*(Fx*ux + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*(uz*uz) + 0.3333333333333333*(ux*ux + uy*uy + (-2. + uz)*uz)) - 
   0.0625*(dirGradC5 - nx*ux - ny*uy - nz*uz)*(2*chem*(uz*uz) - 0.3333333333333333*((-rhoA + rhoB)*(uz*uz) + 2*chem*(ux*ux + uy*uy + (-2 + uz)*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux + uy*uy + (-2 + uz)*uz)));  
                feq6 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-(uz*uz) + 0.3333333333333333*(ux*ux + uy*uy + uz*(2 + uz))) - 
   0.125*(Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(uz*uz) + 0.3333333333333333*(ux*ux + uy*uy + uz*(2. + uz))) - 
   0.0625*(dirGradC6 - nx*ux - ny*uy - nz*uz)*(2*chem*(uz*uz) - 0.3333333333333333*((-rhoA + rhoB)*(uz*uz) + 2*chem*(ux*ux + uy*uy + uz*(2 + uz))) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux + uy*uy + uz*(2 + uz)))); 
                feq7 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((ux + uy)*(ux + uy)) + 0.3333333333333333*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) - 
   0.0625*(Fx*(-1. + ux) + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*(ux*ux) - 2.*ux*uy - 1.*(uy*uy) + 
      0.3333333333333333*(-2.*ux + ux*ux - 2.*uy + uy*uy + uz*uz)) - 
   0.03125*(dirGradC7 - nx*ux - ny*uy - nz*uz)*(2*chem*((ux + uy)*(ux + uy)) + 
      0.3333333333333333*((rhoA - rhoB)*((ux + uy)*(ux + uy)) - 2*chem*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz))); 
                feq8 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((ux + uy)*(ux + uy)) + 0.3333333333333333*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) - 
   0.0625*(Fx + Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(ux*ux) - 2.*ux*uy - 1.*(uy*uy) + 
      0.3333333333333333*(2.*ux + ux*ux + 2.*uy + uy*uy + uz*uz)) - 
   0.03125*(dirGradC8 - nx*ux - ny*uy - nz*uz)*(2*chem*((ux + uy)*(ux + uy)) - 
      0.3333333333333333*(-((rhoA - rhoB)*((ux + uy)*(ux + uy))) + 2*chem*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz))); 
                feq9 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((ux - uy)*(ux - uy)) + 0.3333333333333333*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) - 
   0.0625*(Fy + Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(ux*ux) + 2.*ux*uy - 1.*(uy*uy) + 
      0.3333333333333333*(-2.*ux + ux*ux + 2.*uy + uy*uy + uz*uz)) - 
   0.03125*(dirGradC9 - nx*ux - ny*uy - nz*uz)*(2*chem*((ux - uy)*(ux - uy)) - 
      0.3333333333333333*(-((rhoA - rhoB)*((ux - uy)*(ux - uy))) + 2*chem*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz))); 
                feq10 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((ux - uy)*(ux - uy)) + 0.3333333333333333*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) - 
   0.0625*(Fx*(1 + ux) + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*(ux*ux) + 2.*ux*uy - 1.*(uy*uy) + 
      0.3333333333333333*(2.*ux + ux*ux - 2.*uy + uy*uy + uz*uz)) - 
   0.03125*(dirGradC10 - nx*ux - ny*uy - nz*uz)*(2*chem*((ux - uy)*(ux - uy)) - 
      0.3333333333333333*(-((rhoA - rhoB)*((ux - uy)*(ux - uy))) + 2*chem*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz))); 
                feq11 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((ux + uz)*(ux + uz)) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) - 
   0.0625*(Fx*(-1. + ux) + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*(ux*ux) - 2.*ux*uz - 1.*(uz*uz) + 
      0.3333333333333333*(-2.*ux + ux*ux + uy*uy + (-2. + uz)*uz)) - 
   0.03125*(dirGradC11 - nx*ux - ny*uy - nz*uz)*(2*chem*((ux + uz)*(ux + uz)) + 
      0.3333333333333333*((rhoA - rhoB)*((ux + uz)*(ux + uz)) - 2*chem*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz))); 
                feq12 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((ux + uz)*(ux + uz)) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) - 
   0.0625*(Fx + Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(ux*ux) - 2.*ux*uz - 1.*(uz*uz) + 
      0.3333333333333333*(2.*ux + ux*ux + uy*uy + uz*(2. + uz))) - 
   0.03125*(dirGradC12 - nx*ux - ny*uy - nz*uz)*(2*chem*((ux + uz)*(ux + uz)) - 
      0.3333333333333333*(-((rhoA - rhoB)*((ux + uz)*(ux + uz))) + 2*chem*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))); 
                feq13 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((ux - uz)*(ux - uz)) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) - 
   0.0625*(Fz + Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(ux*ux) + 2.*ux*uz - 1.*(uz*uz) + 
      0.3333333333333333*(-2.*ux + ux*ux + uy*uy + uz*(2. + uz))) - 
   0.03125*(dirGradC13 - nx*ux - ny*uy - nz*uz)*(2*chem*((ux - uz)*(ux - uz)) - 
      0.3333333333333333*(-((rhoA - rhoB)*((ux - uz)*(ux - uz))) + 2*chem*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*(2 + uz)))); 
                feq14 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((ux - uz)*(ux - uz)) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) - 
   0.0625*(Fx*(1 + ux) + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*(ux*ux) + 2.*ux*uz - 1.*(uz*uz) + 
      0.3333333333333333*(2.*ux + ux*ux + uy*uy + (-2. + uz)*uz)) - 
   0.03125*(dirGradC14 - nx*ux - ny*uy - nz*uz)*(2*chem*((ux - uz)*(ux - uz)) - 
      0.3333333333333333*(-((rhoA - rhoB)*((ux - uz)*(ux - uz))) + 2*chem*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz))); 
                feq15 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((uy + uz)*(uy + uz)) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) - 
   0.0625*(Fx*ux + Fy*(-1. + uy) + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*(uy*uy) - 2.*uy*uz - 1.*(uz*uz) + 
      0.3333333333333333*(ux*ux - 2.*uy + uy*uy + (-2. + uz)*uz)) - 
   0.03125*(dirGradC15 - nx*ux - ny*uy - nz*uz)*(2*chem*((uy + uz)*(uy + uz)) + 
      0.3333333333333333*((rhoA - rhoB)*((uy + uz)*(uy + uz)) - 2*chem*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))); 
                feq16 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((uy + uz)*(uy + uz)) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) - 
   0.0625*(Fy + Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(uy*uy) - 2.*uy*uz - 1.*(uz*uz) + 
      0.3333333333333333*(ux*ux + 2.*uy + uy*uy + uz*(2. + uz))) - 
   0.03125*(dirGradC16 - nx*ux - ny*uy - nz*uz)*(2*chem*((uy + uz)*(uy + uz)) - 
      0.3333333333333333*(-((rhoA - rhoB)*((uy + uz)*(uy + uz))) + 2*chem*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))); 
                feq17 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((uy - uz)*(uy - uz)) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) - 
   0.0625*(Fz + Fx*ux + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*(uy*uy) + 2.*uy*uz - 1.*(uz*uz) + 
      0.3333333333333333*(ux*ux - 2.*uy + uy*uy + uz*(2. + uz))) - 
   0.03125*(dirGradC17 - nx*ux - ny*uy - nz*uz)*(2*chem*((uy - uz)*(uy - uz)) - 
      0.3333333333333333*(-((rhoA - rhoB)*((uy - uz)*(uy - uz))) + 2*chem*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))); 
                feq18 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((uy - uz)*(uy - uz)) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) - 
   0.0625*(Fx*ux + Fy*(1 + uy) + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*(uy*uy) + 2.*uy*uz - 1.*(uz*uz) + 
      0.3333333333333333*(ux*ux + 2.*uy + uy*uy + (-2. + uz)*uz)) - 
   0.03125*(dirGradC18 - nx*ux - ny*uy - nz*uz)*(2*chem*((uy - uz)*(uy - uz)) - 
      0.3333333333333333*(-((rhoA - rhoB)*((uy - uz)*(uy - uz))) + 2*chem*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz))); 

            //------------------------------------------------- BCK collison ------------------------------------------------------------//
    		// q=0
    		dist[n] = m0 - (m0-feq0)/tau + 0.25*(2*(Fx*ux + Fy*uy + Fz*uz)*(-0.6666666666666666 + ux*ux + uy*uy + uz*uz) + 
     (mgx*ux + mgy*uy + mgz*uz)*(2*chem*(ux*ux + uy*uy + uz*uz) + 0.3333333333333333*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + uz*uz))));

    		// q = 1
    		dist[nr2] = m1 - (m1-feq1)/tau + 0.125*(2*(Fx*(-1 + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - ux*ux + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*uz)) - 
     (dirGradM1 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*(ux*ux) + 0.3333333333333333*((-rhoA + rhoB)*(ux*ux) + 2*chem*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*uz))));

    		// q=2
    		dist[nr1] = m2 - (m2-feq2)/tau + 0.125*(2*(Fx + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - ux*ux + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*uz)) - 
     (dirGradM2 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*(ux*ux) + 0.3333333333333333*((-rhoA + rhoB)*(ux*ux) + 2*chem*(2*ux + ux*ux + uy*uy + uz*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*uz))));

    		// q = 3
    		dist[nr4] = m3 - (m3-feq3)/tau + 0.125*(2*(Fx*ux + Fy*(-1 + uy) + Fz*uz)*(-0.2222222222222222 - uy*uy + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*uz)) - 
     (dirGradM3 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*(uy*uy) + 0.3333333333333333*((-rhoA + rhoB)*(uy*uy) + 2*chem*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*uz))));

    		// q = 4
    		dist[nr3] = m4 - (m4-feq4)/tau + 0.125*(2*(Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - uy*uy + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*uz)) - 
     (dirGradM4 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*(uy*uy) + 0.3333333333333333*((-rhoA + rhoB)*(uy*uy) + 2*chem*(ux*ux + 2*uy + uy*uy + uz*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*uz))));

    		// q = 5
    		dist[nr6] = m5 - (m5-feq5)/tau + 0.125*(2*(Fx*ux + Fy*uy + Fz*(-1 + uz))*(-0.2222222222222222 - uz*uz + 0.3333333333333333*(ux*ux + uy*uy + (-2 + uz)*uz)) - 
     (dirGradM5 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*(uz*uz) + 0.3333333333333333*((-rhoA + rhoB)*(uz*uz) + 2*chem*(ux*ux + uy*uy + (-2 + uz)*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + (-2 + uz)*uz))));

    		// q = 6
    		dist[nr5] = m6 - (m6-feq6)/tau + 0.125*(2*(Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - uz*uz + 0.3333333333333333*(ux*ux + uy*uy + uz*(2 + uz))) - 
     (dirGradM6 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*(uz*uz) + 0.3333333333333333*((-rhoA + rhoB)*(uz*uz) + 2*chem*(ux*ux + uy*uy + uz*(2 + uz))) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + uz*(2 + uz)))));

    		// q = 7
    		dist[nr8] = m7 - (m7-feq7)/tau + 0.0625*(-2*(Fx*(-1 + ux) + Fy*(-1 + uy) + Fz*uz)*(0.2222222222222222 + (ux + uy)*(ux + uy) - 0.3333333333333333*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) - 
     (dirGradM7 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((ux + uy)*(ux + uy)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((ux + uy)*(ux + uy))) + 2*chem*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz))));

    		// q = 8
    		dist[nr7] = m8 - (m8-feq8)/tau + 0.0625*(2*(Fx + Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - (ux + uy)*(ux + uy) + 0.3333333333333333*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) - 
     (dirGradM8 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((ux + uy)*(ux + uy)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((ux + uy)*(ux + uy))) + 2*chem*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz))));

    		// q = 9
    		dist[nr10] = m9 - (m9-feq9)/tau + 0.0625*(2*(Fy + Fx*(-1 + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - (ux - uy)*(ux - uy) + 0.3333333333333333*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) - 
     (dirGradM9 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((ux - uy)*(ux - uy)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((ux - uy)*(ux - uy))) + 2*chem*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz))));

    		// q = 10
    		dist[nr9] = m10 - (m10-feq10)/tau + 0.0625*(2*(Fx*(1 + ux) + Fy*(-1 + uy) + Fz*uz)*(-0.2222222222222222 - (ux - uy)*(ux - uy) + 0.3333333333333333*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) - 
     (dirGradM10 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((ux - uy)*(ux - uy)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((ux - uy)*(ux - uy))) + 2*chem*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz))));

    		// q = 11
    		dist[nr12] = m11 - (m11-feq11)/tau + 0.0625*(-2*(Fx*(-1 + ux) + Fy*uy + Fz*(-1 + uz))*(0.2222222222222222 + (ux + uz)*(ux + uz) - 0.3333333333333333*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) - 
     (dirGradM11 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((ux + uz)*(ux + uz)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((ux + uz)*(ux + uz))) + 2*chem*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz))));

    		// q = 12
    		dist[nr11] = m12 - (m12-feq12)/tau + 0.0625*(2*(Fx + Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - (ux + uz)*(ux + uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) - 
     (dirGradM12 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((ux + uz)*(ux + uz)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((ux + uz)*(ux + uz))) + 2*chem*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))));

    		// q = 13
    		dist[nr14] = m13 - (m13-feq13)/tau + 0.0625*(2*(Fz + Fx*(-1 + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - (ux - uz)*(ux - uz) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) - 
     (dirGradM13 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((ux - uz)*(ux - uz)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((ux - uz)*(ux - uz))) + 2*chem*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*(2 + uz)))));

    		// q= 14
    		dist[nr13] = m14 - (m14-feq14)/tau + 0.0625*(2*(Fx*(1 + ux) + Fy*uy + Fz*(-1 + uz))*(-0.2222222222222222 - (ux - uz)*(ux - uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) - 
     (dirGradM14 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((ux - uz)*(ux - uz)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((ux - uz)*(ux - uz))) + 2*chem*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz))));

    		// q = 15
    		dist[nr16] = m15 - (m15-feq15)/tau + 0.0625*(-2*(Fx*ux + Fy*(-1 + uy) + Fz*(-1 + uz))*(0.2222222222222222 + (uy + uz)*(uy + uz) - 0.3333333333333333*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) - 
     (dirGradM15 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((uy + uz)*(uy + uz)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((uy + uz)*(uy + uz))) + 2*chem*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))));

    		// q = 16
    		dist[nr15] = m16 - (m16-feq16)/tau + 0.0625*(2*(Fy + Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - (uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) - 
     (dirGradM16 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((uy + uz)*(uy + uz)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((uy + uz)*(uy + uz))) + 2*chem*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))));

    		// q = 17
    		dist[nr18] = m17 - (m17-feq17)/tau + 0.0625*(2*(Fz + Fx*ux + Fy*(-1 + uy) + Fz*uz)*(-0.2222222222222222 - (uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) - 
     (dirGradM17 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((uy - uz)*(uy - uz)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((uy - uz)*(uy - uz))) + 2*chem*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))));

    		// q = 18
    		dist[nr17] = m18 - (m18-feq18)/tau + 0.0625*(2*(Fx*ux + Fy*(1 + uy) + Fz*(-1 + uz))*(-0.2222222222222222 - (uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) - 
     (dirGradM18 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((uy - uz)*(uy - uz)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((uy - uz)*(uy - uz))) + 2*chem*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz))));
            //----------------------------------------------------------------------------------------------------------------------------------------//

            // ----------------------------- compute phase field evolution ----------------------------------------
            //Normalize the Color Gradient
            C = sqrt(nx*nx+ny*ny+nz*nz);
            double ColorMag = C;
            if (C==0.0) ColorMag=1.0;
            nx = nx/ColorMag;
            ny = ny/ColorMag;
            nz = nz/ColorMag;		
            //compute surface tension-related parameter
            //theta = 4.5*M*2.0*(1-phi*phi)/W;
            theta = 4.5*M*2.0*(1-phi_temp*phi_temp)/W;

            //load distributions of phase field
            //q=0
            h0 = hq[n];
            //q=1
            h1 = hq[nr1]; 

            //q=2
            h2 = hq[nr2];  

            //q=3
            h3 = hq[nr3];

            //q=4
            h4 = hq[nr4];

            //q=5
            h5 = hq[nr5];

            //q=6
            h6 = hq[nr6];

            //-------------------------------- BGK collison for phase field ---------------------------------//
            // q = 0
            hq[n] = h0 - (h0 - 0.3333333333333333*phi)/tauM;

            // q = 1
            hq[nr2] = h1 - (h1 - 0.1111111111111111*nx*theta - phi*(0.1111111111111111 + 0.5*ux))/tauM;

            // q = 2
            hq[nr1] = h2 - (h2 + 0.1111111111111111*nx*theta - phi*(0.1111111111111111 - 0.5*ux))/tauM;

            // q = 3
            hq[nr4] = h3 - (h3 - 0.1111111111111111*ny*theta - phi*(0.1111111111111111 + 0.5*uy))/tauM;

            // q = 4
            hq[nr3] = h4 - (h4 + 0.1111111111111111*ny*theta - phi*(0.1111111111111111 - 0.5*uy))/tauM;

            // q = 5
            hq[nr6] = h5 - (h5 - 0.1111111111111111*nz*theta - phi*(0.1111111111111111 + 0.5*uz))/tauM;

            // q = 6
            hq[nr5] = h6 - (h6 + 0.1111111111111111*nz*theta - phi*(0.1111111111111111 - 0.5*uz))/tauM;
            //........................................................................

            //Update velocity on device
    		Vel[0*Np+n] = ux;
    		Vel[1*Np+n] = uy;
    		Vel[2*Np+n] = uz;
            //Update pressure on device
            Pressure[n] = p;
            //Update chemical potential on device
            mu_phi[n] = chem;
            //Update color gradient on device
    		ColorGrad[0*Np+n] = nx;
    		ColorGrad[1*Np+n] = ny;
    		ColorGrad[2*Np+n] = nz;
    	}
    }
}

__global__ void dvc_ScaLBL_D3Q19_AAeven_FreeLeeModel_Combined(int *Map, double *dist, double *hq, double *Den,	double *Phi, double *mu_phi, double *Vel, double *Pressure, double *ColorGrad,
        double rhoA, double rhoB, double tauA, double tauB, double tauM, double kappa, double beta, double W, double Fx, double Fy, double Fz, 
        int strideY, int strideZ, int start, int finish, int Np){

	int n,nn,nn2x,ijk;
    double ux,uy,uz;//fluid velocity 
    double p;//pressure
    double chem;//chemical potential
    double phi; //phase field
    double rho0;//fluid density
	// distribution functions
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
	double m0,m3,m5,m7;
	//double mm1,mm2,mm4,mm6,mm8,mm9,mm10,mm11,mm12,mm13,mm14,mm15,mm16,mm17,mm18;
	//double mm3,mm5,mm7;
    double feq0,feq1,feq2,feq3,feq4,feq5,feq6,feq7,feq8,feq9,feq10,feq11,feq12,feq13,feq14,feq15,feq16,feq17,feq18;
    double nx,ny,nz;//normal color gradient
    double mgx,mgy,mgz;//mixed gradient reaching secondary neighbor
    double dirGradC1,dirGradC2,dirGradC3,dirGradC4,dirGradC5,dirGradC6,dirGradC7,dirGradC8,dirGradC9,dirGradC10,dirGradC11,dirGradC12;
    double dirGradC13,dirGradC14,dirGradC15,dirGradC16,dirGradC17,dirGradC18;
    double dirGradM1,dirGradM2,dirGradM3,dirGradM4,dirGradM5,dirGradM6,dirGradM7,dirGradM8,dirGradM9,dirGradM10,dirGradM11,dirGradM12;
    double dirGradM13,dirGradM14,dirGradM15,dirGradM16,dirGradM17,dirGradM18;

    double h0,h1,h2,h3,h4,h5,h6;//distributions for LB phase field
	double tau;//position dependent LB relaxation time for fluid
    double C,theta;
    double M = 2.0/9.0*(tauM-0.5);//diffusivity (or mobility) for the phase field D3Q7
    double phi_temp;

    int S = Np/NBLOCKS/NTHREADS + 1;
    for (int s=0; s<S; s++){
    	//........Get 1-D index for this thread....................
    	n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

    	if ( n<finish ){
    		rho0 = Den[n];//load density

    		// Get the 1D index based on regular data layout
    		ijk = Map[n];
            phi = Phi[ijk];// load phase field
            phi_temp = phi;
            if (phi>1.f) phi_temp=1.0;
            if (phi<-1.f) phi_temp=-1.0;

    		// local relaxation time
    		tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);

    		//					COMPUTE THE COLOR GRADIENT
    		//........................................................................
    		//.................Read Phase Indicator Values............................
    		//........................................................................
    		nn = ijk-1;							// neighbor index (get convention)
    		m2 = Phi[nn];						// get neighbor for phi - 1
    		//........................................................................
    		nn = ijk+1;							// neighbor index (get convention)
    		m1 = Phi[nn];						// get neighbor for phi - 2
            dirGradC1 = 0.5*(m1-m2);
            dirGradC2 = 0.5*(m2-m1);
    		//........................................................................
    		nn = ijk-strideY;							// neighbor index (get convention)
    		m4 = Phi[nn];					// get neighbor for phi - 3
    		//........................................................................
    		nn = ijk+strideY;							// neighbor index (get convention)
    		m3 = Phi[nn];					// get neighbor for phi - 4
            dirGradC3 = 0.5*(m3-m4);
            dirGradC4 = 0.5*(m4-m3);
    		//........................................................................
    		nn = ijk-strideZ;						// neighbor index (get convention)
    		m6 = Phi[nn];					// get neighbor for phi - 5
    		//........................................................................
    		nn = ijk+strideZ;						// neighbor index (get convention)
    		m5 = Phi[nn];					// get neighbor for phi - 6
            dirGradC5 = 0.5*(m5-m6);
            dirGradC6 = 0.5*(m6-m5);
    		//........................................................................
    		nn = ijk-strideY-1;						// neighbor index (get convention)
    		m8 = Phi[nn];					// get neighbor for phi - 7
    		//........................................................................
    		nn = ijk+strideY+1;						// neighbor index (get convention)
    		m7 = Phi[nn];					// get neighbor for phi - 8
            dirGradC7 = 0.5*(m7-m8);
            dirGradC8 = 0.5*(m8-m7);
    		//........................................................................
    		nn = ijk+strideY-1;						// neighbor index (get convention)
    		m10 = Phi[nn];					// get neighbor for phi - 9
    		//........................................................................
    		nn = ijk-strideY+1;						// neighbor index (get convention)
    		m9 = Phi[nn];					// get neighbor for phi - 10
            dirGradC9  = 0.5*(m9-m10);
            dirGradC10 = 0.5*(m10-m9);
    		//........................................................................
    		nn = ijk-strideZ-1;						// neighbor index (get convention)
    		m12 = Phi[nn];					// get neighbor for phi - 11
    		//........................................................................
    		nn = ijk+strideZ+1;						// neighbor index (get convention)
    		m11 = Phi[nn];					// get neighbor for phi - 12
            dirGradC11 = 0.5*(m11-m12);
            dirGradC12 = 0.5*(m12-m11);
    		//........................................................................
    		nn = ijk+strideZ-1;						// neighbor index (get convention)
    		m14 = Phi[nn];					// get neighbor for phi - 13
    		//........................................................................
    		nn = ijk-strideZ+1;						// neighbor index (get convention)
    		m13 = Phi[nn];					// get neighbor for phi - 14
            dirGradC13 = 0.5*(m13-m14);
            dirGradC14 = 0.5*(m14-m13);
    		//........................................................................
    		nn = ijk-strideZ-strideY;					// neighbor index (get convention)
    		m16 = Phi[nn];					// get neighbor for phi - 15
    		//........................................................................
    		nn = ijk+strideZ+strideY;					// neighbor index (get convention)
    		m15 = Phi[nn];					// get neighbor for phi - 16
            dirGradC15 = 0.5*(m15-m16);
            dirGradC16 = 0.5*(m16-m15);
    		//........................................................................
    		nn = ijk+strideZ-strideY;					// neighbor index (get convention)
    		m18 = Phi[nn];					// get neighbor for phi - 17
    		//........................................................................
    		nn = ijk-strideZ+strideY;					// neighbor index (get convention)
    		m17 = Phi[nn];					// get neighbor for phi - 18
            dirGradC17 = 0.5*(m17-m18);
            dirGradC18 = 0.5*(m18-m17);

            // compute mixed difference (Eq.30, A.Fukhari et al. JCP 315(2016) 434-457)
    		//........................................................................
    		nn2x = ijk+2;							// neighbor index (get convention)
    		dirGradM1 = Phi[nn2x];						// get neighbor for phi - 1
            dirGradM1 = 0.25*(-dirGradM1+5.0*m1-3.0*phi-m2);
    		//........................................................................
    		nn2x = ijk-2;							// neighbor index (get convention)
    		dirGradM2 = Phi[nn2x];						// get neighbor for phi - 2
            dirGradM2 = 0.25*(-dirGradM2+5.0*m2-3.0*phi-m1);
    		//........................................................................
    		nn2x = ijk+strideY*2;							// neighbor index (get convention)
    		dirGradM3 = Phi[nn2x];					// get neighbor for phi - 3
            dirGradM3 = 0.25*(-dirGradM3+5.0*m3-3.0*phi-m4);
    		//........................................................................
    		nn2x = ijk-strideY*2;							// neighbor index (get convention)
    		dirGradM4 = Phi[nn2x];					// get neighbor for phi - 4
            dirGradM4 = 0.25*(-dirGradM4+5.0*m4-3.0*phi-m3);
    		//........................................................................
    		nn2x = ijk+strideZ*2;						// neighbor index (get convention)
    		dirGradM5 = Phi[nn2x];					// get neighbor for phi - 5
            dirGradM5 = 0.25*(-dirGradM5+5.0*m5-3.0*phi-m6);
    		//........................................................................
    		nn2x = ijk-strideZ*2;						// neighbor index (get convention)
    		dirGradM6 = Phi[nn2x];					// get neighbor for phi - 6
            dirGradM6 = 0.25*(-dirGradM6+5.0*m6-3.0*phi-m5);
    		//........................................................................
    		nn2x = ijk+strideY*2+2;						// neighbor index (get convention)
    		dirGradM7 = Phi[nn2x];					// get neighbor for phi - 7
            dirGradM7 = 0.25*(-dirGradM7+5.0*m7-3.0*phi-m8);
    		//........................................................................
    		nn2x = ijk-strideY*2-2;						// neighbor index (get convention)
    		dirGradM8 = Phi[nn2x];					// get neighbor for phi - 8
            dirGradM8 = 0.25*(-dirGradM8+5.0*m8-3.0*phi-m7);
    		//........................................................................
    		nn2x = ijk-strideY*2+2;						// neighbor index (get convention)
    		dirGradM9 = Phi[nn2x];					// get neighbor for phi - 9
            dirGradM9 = 0.25*(-dirGradM9+5.0*m9-3.0*phi-m10);
    		//........................................................................
    		nn2x = ijk+strideY*2-2;						// neighbor index (get convention)
    		dirGradM10 = Phi[nn2x];					// get neighbor for phi - 10
            dirGradM10 = 0.25*(-dirGradM10+5.0*m10-3.0*phi-m9);
    		//........................................................................
    		nn2x = ijk+strideZ*2+2;						// neighbor index (get convention)
    		dirGradM11 = Phi[nn2x];					// get neighbor for phi - 11
            dirGradM11 = 0.25*(-dirGradM11+5.0*m11-3.0*phi-m12);
    		//........................................................................
    		nn2x = ijk-strideZ*2-2;						// neighbor index (get convention)
    		dirGradM12 = Phi[nn2x];					// get neighbor for phi - 12
            dirGradM12 = 0.25*(-dirGradM12+5.0*m12-3.0*phi-m11);
    		//........................................................................
    		nn2x = ijk-strideZ*2+2;						// neighbor index (get convention)
    		dirGradM13 = Phi[nn2x];					// get neighbor for phi - 13
            dirGradM13 = 0.25*(-dirGradM13+5.0*m13-3.0*phi-m14);
    		//........................................................................
    		nn2x = ijk+strideZ*2-2;						// neighbor index (get convention)
    		dirGradM14 = Phi[nn2x];					// get neighbor for phi - 14
            dirGradM14 = 0.25*(-dirGradM14+5.0*m14-3.0*phi-m13);
    		//........................................................................
    		nn2x = ijk+strideZ*2+strideY*2;					// neighbor index (get convention)
    		dirGradM15 = Phi[nn2x];					// get neighbor for phi - 15
            dirGradM15 = 0.25*(-dirGradM15+5.0*m15-3.0*phi-m16);
    		//........................................................................
    		nn2x = ijk-strideZ*2-strideY*2;					// neighbor index (get convention)
    		dirGradM16 = Phi[nn2x];					// get neighbor for phi - 16
            dirGradM16 = 0.25*(-dirGradM16+5.0*m16-3.0*phi-m15);
    		//........................................................................
    		nn2x = ijk-strideZ*2+strideY*2;					// neighbor index (get convention)
    		dirGradM17 = Phi[nn2x];					// get neighbor for phi - 17
            dirGradM17 = 0.25*(-dirGradM17+5.0*m17-3.0*phi-m18);
    		//........................................................................
    		nn2x = ijk+strideZ*2-strideY*2;					// neighbor index (get convention)
    		dirGradM18 = Phi[nn2x];					// get neighbor for phi - 18
            dirGradM18 = 0.25*(-dirGradM18+5.0*m18-3.0*phi-m17);


    		//............Compute the Color Gradient...................................
    		nx = 3.0*1.0/18.0*(dirGradC1-dirGradC2+0.5*(dirGradC7-dirGradC8+dirGradC9-dirGradC10+dirGradC11-dirGradC12+dirGradC13-dirGradC14));
    		ny = 3.0*1.0/18.0*(dirGradC3-dirGradC4+0.5*(dirGradC7-dirGradC8-dirGradC9+dirGradC10+dirGradC15-dirGradC16+dirGradC17-dirGradC18));
    		nz = 3.0*1.0/18.0*(dirGradC5-dirGradC6+0.5*(dirGradC11-dirGradC12-dirGradC13+dirGradC14+dirGradC15-dirGradC16-dirGradC17+dirGradC18));
    		//............Compute the Chemical Potential...............................
            //chem = 2.0*3.0/18.0*(m1+m2+m3+m4+m5+m6-6*phi+0.5*(m7+m8+m9+m10+m11+m12+m13+m14+m15+m16+m17+m18-12*phi));//intermediate var, i.e. the laplacian
            //chem = 4.0*beta*phi*(phi+1.0)*(phi-1.0)-kappa*chem;
            chem = 2.0*3.0/18.0*(m1+m2+m3+m4+m5+m6-6*phi_temp+0.5*(m7+m8+m9+m10+m11+m12+m13+m14+m15+m16+m17+m18-12*phi_temp));//intermediate var, i.e. the laplacian
            chem = 4.0*beta*phi_temp*(phi_temp+1.0)*(phi_temp-1.0)-kappa*chem;
    		//............Compute the Mixed Gradient...................................
    		mgx = 3.0*1.0/18.0*(dirGradM1-dirGradM2+0.5*(dirGradM7-dirGradM8+dirGradM9-dirGradM10+dirGradM11-dirGradM12+dirGradM13-dirGradM14));
    		mgy = 3.0*1.0/18.0*(dirGradM3-dirGradM4+0.5*(dirGradM7-dirGradM8-dirGradM9+dirGradM10+dirGradM15-dirGradM16+dirGradM17-dirGradM18));
    		mgz = 3.0*1.0/18.0*(dirGradM5-dirGradM6+0.5*(dirGradM11-dirGradM12-dirGradM13+dirGradM14+dirGradM15-dirGradM16-dirGradM17+dirGradM18));

            //de-noise color gradient and mixed gradient
            C = sqrt(nx*nx+ny*ny+nz*nz);
            if (C<1.0e-12) nx=ny=nz=0.0;
            double mg_mag = sqrt(mgx*mgx+mgy*mgy+mgz*mgz);
            if (mg_mag<1.0e-12) mgx=mgy=mgz=0.0;
            //maybe you can also de-noise chemical potential ? within the bulk phase chem should be ZERO
            if (fabs(chem)<1.0e-12) chem=0.0;
    		
    		// q=0
    		m0 = dist[n];
    		// q=1
    		m1 = dist[2*Np+n]; 

            // q=2
    		m2 = dist[1*Np+n];  

    		// q=3
    		m3 = dist[4*Np+n];

    		// q = 4
    		m4 = dist[3*Np+n];

    		// q=5
    		m5 = dist[6*Np+n];

    		// q = 6
    		m6 = dist[5*Np+n];
    		
    		// q=7
    		m7 = dist[8*Np+n];

    		// q = 8
    		m8 = dist[7*Np+n];

    		// q=9
    		m9 = dist[10*Np+n];

    		// q = 10
    		m10 = dist[9*Np+n];

    		// q=11
    		m11 = dist[12*Np+n];

    		// q=12
    		m12 = dist[11*Np+n];

    		// q=13
    		m13 = dist[14*Np+n];

    		// q=14
    		m14 = dist[13*Np+n];

    		// q=15
    		m15 = dist[16*Np+n];

    		// q=16
    		m16 = dist[15*Np+n];

    		// q=17
    		m17 = dist[18*Np+n];

    		// q=18
    		m18 = dist[17*Np+n];

            //compute fluid velocity
            ux = 3.0/rho0*(m1-m2+m7-m8+m9-m10+m11-m12+m13-m14+0.5*(chem*nx+Fx)/3.0);
            uy = 3.0/rho0*(m3-m4+m7-m8-m9+m10+m15-m16+m17-m18+0.5*(chem*ny+Fy)/3.0);
            uz = 3.0/rho0*(m5-m6+m11-m12-m13+m14+m15-m16-m17+m18+0.5*(chem*nz+Fz)/3.0);

            //compute pressure
            p = (m0+m2+m1+m4+m3+m6+m5+m8+m7+m10+m9+m12+m11+m14+m13+m16+m15+m18+m17)
                      +0.5*(rhoA-rhoB)/2.0/3.0*(ux*nx+uy*ny+uz*nz);

            //compute equilibrium distributions
                feq0 = 0.3333333333333333*p - 0.25*(Fx*ux + Fy*uy + Fz*uz)*(-0.6666666666666666 + ux*ux + uy*uy + uz*uz) - 0.16666666666666666*rho0*(ux*ux + uy*uy + uz*uz) - 
   0.5*(-(nx*ux) - ny*uy - nz*uz)*(-0.08333333333333333*(rhoA - rhoB)*(ux*ux + uy*uy + uz*uz) + chem*(0.3333333333333333 - 0.5*(ux*ux + uy*uy + uz*uz))); 
                feq1 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-(ux*ux) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*uz)) - 
   0.125*(Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(ux*ux) + 0.3333333333333333*(-2.*ux + ux*ux + uy*uy + uz*uz)) - 
   0.0625*(dirGradC1 - nx*ux - ny*uy - nz*uz)*(2*chem*(ux*ux) - 0.3333333333333333*((-rhoA + rhoB)*(ux*ux) + 2*chem*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*uz))); 
                feq2 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-(ux*ux) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*uz)) - 
   0.125*(Fx + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(ux*ux) + 0.3333333333333333*(2.*ux + ux*ux + uy*uy + uz*uz)) - 
   0.0625*(dirGradC2 - nx*ux - ny*uy - nz*uz)*(2*chem*(ux*ux) - 0.3333333333333333*((-rhoA + rhoB)*(ux*ux) + 2*chem*(2*ux + ux*ux + uy*uy + uz*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*uz))); 
                feq3 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-(uy*uy) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*uz)) - 
   0.125*(Fx*ux + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*(uy*uy) + 0.3333333333333333*(ux*ux - 2.*uy + uy*uy + uz*uz)) - 
   0.0625*(dirGradC3 - nx*ux - ny*uy - nz*uz)*(2*chem*(uy*uy) - 0.3333333333333333*((-rhoA + rhoB)*(uy*uy) + 2*chem*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*uz)));  
                feq4 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-(uy*uy) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*uz)) - 
   0.125*(Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(uy*uy) + 0.3333333333333333*(ux*ux + 2.*uy + uy*uy + uz*uz)) - 
   0.0625*(dirGradC4 - nx*ux - ny*uy - nz*uz)*(2*chem*(uy*uy) - 0.3333333333333333*((-rhoA + rhoB)*(uy*uy) + 2*chem*(ux*ux + 2*uy + uy*uy + uz*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*uz)));  
                feq5 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-(uz*uz) + 0.3333333333333333*(ux*ux + uy*uy + (-2 + uz)*uz)) - 
   0.125*(Fx*ux + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*(uz*uz) + 0.3333333333333333*(ux*ux + uy*uy + (-2. + uz)*uz)) - 
   0.0625*(dirGradC5 - nx*ux - ny*uy - nz*uz)*(2*chem*(uz*uz) - 0.3333333333333333*((-rhoA + rhoB)*(uz*uz) + 2*chem*(ux*ux + uy*uy + (-2 + uz)*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux + uy*uy + (-2 + uz)*uz)));  
                feq6 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-(uz*uz) + 0.3333333333333333*(ux*ux + uy*uy + uz*(2 + uz))) - 
   0.125*(Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(uz*uz) + 0.3333333333333333*(ux*ux + uy*uy + uz*(2. + uz))) - 
   0.0625*(dirGradC6 - nx*ux - ny*uy - nz*uz)*(2*chem*(uz*uz) - 0.3333333333333333*((-rhoA + rhoB)*(uz*uz) + 2*chem*(ux*ux + uy*uy + uz*(2 + uz))) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux + uy*uy + uz*(2 + uz)))); 
                feq7 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((ux + uy)*(ux + uy)) + 0.3333333333333333*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) - 
   0.0625*(Fx*(-1. + ux) + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*(ux*ux) - 2.*ux*uy - 1.*(uy*uy) + 
      0.3333333333333333*(-2.*ux + ux*ux - 2.*uy + uy*uy + uz*uz)) - 
   0.03125*(dirGradC7 - nx*ux - ny*uy - nz*uz)*(2*chem*((ux + uy)*(ux + uy)) + 
      0.3333333333333333*((rhoA - rhoB)*((ux + uy)*(ux + uy)) - 2*chem*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz))); 
                feq8 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((ux + uy)*(ux + uy)) + 0.3333333333333333*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) - 
   0.0625*(Fx + Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(ux*ux) - 2.*ux*uy - 1.*(uy*uy) + 
      0.3333333333333333*(2.*ux + ux*ux + 2.*uy + uy*uy + uz*uz)) - 
   0.03125*(dirGradC8 - nx*ux - ny*uy - nz*uz)*(2*chem*((ux + uy)*(ux + uy)) - 
      0.3333333333333333*(-((rhoA - rhoB)*((ux + uy)*(ux + uy))) + 2*chem*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz))); 
                feq9 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((ux - uy)*(ux - uy)) + 0.3333333333333333*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) - 
   0.0625*(Fy + Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(ux*ux) + 2.*ux*uy - 1.*(uy*uy) + 
      0.3333333333333333*(-2.*ux + ux*ux + 2.*uy + uy*uy + uz*uz)) - 
   0.03125*(dirGradC9 - nx*ux - ny*uy - nz*uz)*(2*chem*((ux - uy)*(ux - uy)) - 
      0.3333333333333333*(-((rhoA - rhoB)*((ux - uy)*(ux - uy))) + 2*chem*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz))); 
                feq10 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((ux - uy)*(ux - uy)) + 0.3333333333333333*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) - 
   0.0625*(Fx*(1 + ux) + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*(ux*ux) + 2.*ux*uy - 1.*(uy*uy) + 
      0.3333333333333333*(2.*ux + ux*ux - 2.*uy + uy*uy + uz*uz)) - 
   0.03125*(dirGradC10 - nx*ux - ny*uy - nz*uz)*(2*chem*((ux - uy)*(ux - uy)) - 
      0.3333333333333333*(-((rhoA - rhoB)*((ux - uy)*(ux - uy))) + 2*chem*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz))); 
                feq11 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((ux + uz)*(ux + uz)) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) - 
   0.0625*(Fx*(-1. + ux) + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*(ux*ux) - 2.*ux*uz - 1.*(uz*uz) + 
      0.3333333333333333*(-2.*ux + ux*ux + uy*uy + (-2. + uz)*uz)) - 
   0.03125*(dirGradC11 - nx*ux - ny*uy - nz*uz)*(2*chem*((ux + uz)*(ux + uz)) + 
      0.3333333333333333*((rhoA - rhoB)*((ux + uz)*(ux + uz)) - 2*chem*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz))); 
                feq12 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((ux + uz)*(ux + uz)) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) - 
   0.0625*(Fx + Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(ux*ux) - 2.*ux*uz - 1.*(uz*uz) + 
      0.3333333333333333*(2.*ux + ux*ux + uy*uy + uz*(2. + uz))) - 
   0.03125*(dirGradC12 - nx*ux - ny*uy - nz*uz)*(2*chem*((ux + uz)*(ux + uz)) - 
      0.3333333333333333*(-((rhoA - rhoB)*((ux + uz)*(ux + uz))) + 2*chem*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))); 
                feq13 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((ux - uz)*(ux - uz)) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) - 
   0.0625*(Fz + Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(ux*ux) + 2.*ux*uz - 1.*(uz*uz) + 
      0.3333333333333333*(-2.*ux + ux*ux + uy*uy + uz*(2. + uz))) - 
   0.03125*(dirGradC13 - nx*ux - ny*uy - nz*uz)*(2*chem*((ux - uz)*(ux - uz)) - 
      0.3333333333333333*(-((rhoA - rhoB)*((ux - uz)*(ux - uz))) + 2*chem*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*(2 + uz)))); 
                feq14 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((ux - uz)*(ux - uz)) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) - 
   0.0625*(Fx*(1 + ux) + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*(ux*ux) + 2.*ux*uz - 1.*(uz*uz) + 
      0.3333333333333333*(2.*ux + ux*ux + uy*uy + (-2. + uz)*uz)) - 
   0.03125*(dirGradC14 - nx*ux - ny*uy - nz*uz)*(2*chem*((ux - uz)*(ux - uz)) - 
      0.3333333333333333*(-((rhoA - rhoB)*((ux - uz)*(ux - uz))) + 2*chem*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz))); 
                feq15 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((uy + uz)*(uy + uz)) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) - 
   0.0625*(Fx*ux + Fy*(-1. + uy) + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*(uy*uy) - 2.*uy*uz - 1.*(uz*uz) + 
      0.3333333333333333*(ux*ux - 2.*uy + uy*uy + (-2. + uz)*uz)) - 
   0.03125*(dirGradC15 - nx*ux - ny*uy - nz*uz)*(2*chem*((uy + uz)*(uy + uz)) + 
      0.3333333333333333*((rhoA - rhoB)*((uy + uz)*(uy + uz)) - 2*chem*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))); 
                feq16 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((uy + uz)*(uy + uz)) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) - 
   0.0625*(Fy + Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(uy*uy) - 2.*uy*uz - 1.*(uz*uz) + 
      0.3333333333333333*(ux*ux + 2.*uy + uy*uy + uz*(2. + uz))) - 
   0.03125*(dirGradC16 - nx*ux - ny*uy - nz*uz)*(2*chem*((uy + uz)*(uy + uz)) - 
      0.3333333333333333*(-((rhoA - rhoB)*((uy + uz)*(uy + uz))) + 2*chem*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))); 
                feq17 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((uy - uz)*(uy - uz)) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) - 
   0.0625*(Fz + Fx*ux + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*(uy*uy) + 2.*uy*uz - 1.*(uz*uz) + 
      0.3333333333333333*(ux*ux - 2.*uy + uy*uy + uz*(2. + uz))) - 
   0.03125*(dirGradC17 - nx*ux - ny*uy - nz*uz)*(2*chem*((uy - uz)*(uy - uz)) - 
      0.3333333333333333*(-((rhoA - rhoB)*((uy - uz)*(uy - uz))) + 2*chem*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))); 
                feq18 = 0.027777777777777776*p - 0.041666666666666664*rho0*(-((uy - uz)*(uy - uz)) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) - 
   0.0625*(Fx*ux + Fy*(1 + uy) + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*(uy*uy) + 2.*uy*uz - 1.*(uz*uz) + 
      0.3333333333333333*(ux*ux + 2.*uy + uy*uy + (-2. + uz)*uz)) - 
   0.03125*(dirGradC18 - nx*ux - ny*uy - nz*uz)*(2*chem*((uy - uz)*(uy - uz)) - 
      0.3333333333333333*(-((rhoA - rhoB)*((uy - uz)*(uy - uz))) + 2*chem*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 
      0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz))); 

            //------------------------------------------------- BCK collison ------------------------------------------------------------//
    		// q=0
    		dist[n] = m0 - (m0-feq0)/tau + 0.25*(2*(Fx*ux + Fy*uy + Fz*uz)*(-0.6666666666666666 + ux*ux + uy*uy + uz*uz) + 
     (mgx*ux + mgy*uy + mgz*uz)*(2*chem*(ux*ux + uy*uy + uz*uz) + 0.3333333333333333*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + uz*uz))));

    		// q = 1
    		dist[1*Np+n] = m1 - (m1-feq1)/tau + 0.125*(2*(Fx*(-1 + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - ux*ux + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*uz)) - 
     (dirGradM1 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*(ux*ux) + 0.3333333333333333*((-rhoA + rhoB)*(ux*ux) + 2*chem*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*uz))));

    		// q=2
    		dist[2*Np+n] = m2 - (m2-feq2)/tau + 0.125*(2*(Fx + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - ux*ux + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*uz)) - 
     (dirGradM2 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*(ux*ux) + 0.3333333333333333*((-rhoA + rhoB)*(ux*ux) + 2*chem*(2*ux + ux*ux + uy*uy + uz*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*uz))));

    		// q = 3
    		dist[3*Np+n] = m3 - (m3-feq3)/tau + 0.125*(2*(Fx*ux + Fy*(-1 + uy) + Fz*uz)*(-0.2222222222222222 - uy*uy + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*uz)) - 
     (dirGradM3 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*(uy*uy) + 0.3333333333333333*((-rhoA + rhoB)*(uy*uy) + 2*chem*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*uz))));

    		// q = 4
    		dist[4*Np+n] = m4 - (m4-feq4)/tau + 0.125*(2*(Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - uy*uy + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*uz)) - 
     (dirGradM4 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*(uy*uy) + 0.3333333333333333*((-rhoA + rhoB)*(uy*uy) + 2*chem*(ux*ux + 2*uy + uy*uy + uz*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*uz))));

    		// q = 5
    		dist[5*Np+n] = m5 - (m5-feq5)/tau + 0.125*(2*(Fx*ux + Fy*uy + Fz*(-1 + uz))*(-0.2222222222222222 - uz*uz + 0.3333333333333333*(ux*ux + uy*uy + (-2 + uz)*uz)) - 
     (dirGradM5 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*(uz*uz) + 0.3333333333333333*((-rhoA + rhoB)*(uz*uz) + 2*chem*(ux*ux + uy*uy + (-2 + uz)*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + (-2 + uz)*uz))));

    		// q = 6
    		dist[6*Np+n] = m6 - (m6-feq6)/tau + 0.125*(2*(Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - uz*uz + 0.3333333333333333*(ux*ux + uy*uy + uz*(2 + uz))) - 
     (dirGradM6 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*(uz*uz) + 0.3333333333333333*((-rhoA + rhoB)*(uz*uz) + 2*chem*(ux*ux + uy*uy + uz*(2 + uz))) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + uz*(2 + uz)))));

    		// q = 7
    		dist[7*Np+n] = m7 - (m7-feq7)/tau + 0.0625*(-2*(Fx*(-1 + ux) + Fy*(-1 + uy) + Fz*uz)*(0.2222222222222222 + (ux + uy)*(ux + uy) - 0.3333333333333333*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) - 
     (dirGradM7 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((ux + uy)*(ux + uy)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((ux + uy)*(ux + uy))) + 2*chem*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz))));

    		// q = 8
    		dist[8*Np+n] = m8 - (m8-feq8)/tau + 0.0625*(2*(Fx + Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - (ux + uy)*(ux + uy) + 0.3333333333333333*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) - 
     (dirGradM8 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((ux + uy)*(ux + uy)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((ux + uy)*(ux + uy))) + 2*chem*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz))));

    		// q = 9
    		dist[9*Np+n] = m9 - (m9-feq9)/tau + 0.0625*(2*(Fy + Fx*(-1 + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - (ux - uy)*(ux - uy) + 0.3333333333333333*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) - 
     (dirGradM9 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((ux - uy)*(ux - uy)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((ux - uy)*(ux - uy))) + 2*chem*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz))));

    		// q = 10
    		dist[10*Np+n] = m10 - (m10-feq10)/tau + 0.0625*(2*(Fx*(1 + ux) + Fy*(-1 + uy) + Fz*uz)*(-0.2222222222222222 - (ux - uy)*(ux - uy) + 0.3333333333333333*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) - 
     (dirGradM10 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((ux - uy)*(ux - uy)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((ux - uy)*(ux - uy))) + 2*chem*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz))));

    		// q = 11
    		dist[11*Np+n] = m11 - (m11-feq11)/tau + 0.0625*(-2*(Fx*(-1 + ux) + Fy*uy + Fz*(-1 + uz))*(0.2222222222222222 + (ux + uz)*(ux + uz) - 0.3333333333333333*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) - 
     (dirGradM11 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((ux + uz)*(ux + uz)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((ux + uz)*(ux + uz))) + 2*chem*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz))));

    		// q = 12
    		dist[12*Np+n] = m12 - (m12-feq12)/tau + 0.0625*(2*(Fx + Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - (ux + uz)*(ux + uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) - 
     (dirGradM12 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((ux + uz)*(ux + uz)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((ux + uz)*(ux + uz))) + 2*chem*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))));

    		// q = 13
    		dist[13*Np+n] = m13 - (m13-feq13)/tau + 0.0625*(2*(Fz + Fx*(-1 + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - (ux - uz)*(ux - uz) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) - 
     (dirGradM13 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((ux - uz)*(ux - uz)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((ux - uz)*(ux - uz))) + 2*chem*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*(2 + uz)))));

    		// q= 14
    		dist[14*Np+n] = m14 - (m14-feq14)/tau + 0.0625*(2*(Fx*(1 + ux) + Fy*uy + Fz*(-1 + uz))*(-0.2222222222222222 - (ux - uz)*(ux - uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) - 
     (dirGradM14 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((ux - uz)*(ux - uz)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((ux - uz)*(ux - uz))) + 2*chem*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz))));

    		// q = 15
    		dist[15*Np+n] = m15 - (m15-feq15)/tau + 0.0625*(-2*(Fx*ux + Fy*(-1 + uy) + Fz*(-1 + uz))*(0.2222222222222222 + (uy + uz)*(uy + uz) - 0.3333333333333333*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) - 
     (dirGradM15 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((uy + uz)*(uy + uz)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((uy + uz)*(uy + uz))) + 2*chem*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))));

    		// q = 16
    		dist[16*Np+n] = m16 - (m16-feq16)/tau + 0.0625*(2*(Fy + Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - (uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) - 
     (dirGradM16 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((uy + uz)*(uy + uz)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((uy + uz)*(uy + uz))) + 2*chem*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))));

    		// q = 17
    		dist[17*Np+n] = m17 - (m17-feq17)/tau + 0.0625*(2*(Fz + Fx*ux + Fy*(-1 + uy) + Fz*uz)*(-0.2222222222222222 - (uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) - 
     (dirGradM17 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((uy - uz)*(uy - uz)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((uy - uz)*(uy - uz))) + 2*chem*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))));

    		// q = 18
    		dist[18*Np+n] = m18 - (m18-feq18)/tau + 0.0625*(2*(Fx*ux + Fy*(1 + uy) + Fz*(-1 + uz))*(-0.2222222222222222 - (uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) - 
     (dirGradM18 - mgx*ux - mgy*uy - mgz*uz)*(-2*chem*((uy - uz)*(uy - uz)) + 
        0.3333333333333333*(-((rhoA - rhoB)*((uy - uz)*(uy - uz))) + 2*chem*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 
        0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz))));
            //----------------------------------------------------------------------------------------------------------------------------------------//

            // ----------------------------- compute phase field evolution ----------------------------------------
            //Normalize the Color Gradient
            C = sqrt(nx*nx+ny*ny+nz*nz);
            double ColorMag = C;
            if (C==0.0) ColorMag=1.0;
            nx = nx/ColorMag;
            ny = ny/ColorMag;
            nz = nz/ColorMag;		
            //compute surface tension-related parameter
            //theta = 4.5*M*2.0*(1-phi*phi)/W;
            theta = 4.5*M*2.0*(1-phi_temp*phi_temp)/W;

            //load distributions of phase field
            //q=0
            h0 = hq[n];
            //q=1
            h1 = hq[2*Np+n]; 

            //q=2
            h2 = hq[1*Np+n];  

            //q=3
            h3 = hq[4*Np+n];

            //q=4
            h4 = hq[3*Np+n];

            //q=5
            h5 = hq[6*Np+n];

            //q=6
            h6 = hq[5*Np+n];

            //-------------------------------- BGK collison for phase field ---------------------------------//
            // q = 0
            hq[n] = h0 - (h0 - 0.3333333333333333*phi)/tauM;

            // q = 1
            hq[1*Np+n] = h1 - (h1 - 0.1111111111111111*nx*theta - phi*(0.1111111111111111 + 0.5*ux))/tauM;

            // q = 2
            hq[2*Np+n] = h2 - (h2 + 0.1111111111111111*nx*theta - phi*(0.1111111111111111 - 0.5*ux))/tauM;

            // q = 3
            hq[3*Np+n] = h3 - (h3 - 0.1111111111111111*ny*theta - phi*(0.1111111111111111 + 0.5*uy))/tauM;

            // q = 4
            hq[4*Np+n] = h4 - (h4 + 0.1111111111111111*ny*theta - phi*(0.1111111111111111 - 0.5*uy))/tauM;

            // q = 5
            hq[5*Np+n] = h5 - (h5 - 0.1111111111111111*nz*theta - phi*(0.1111111111111111 + 0.5*uz))/tauM;

            // q = 6
            hq[6*Np+n] = h6 - (h6 + 0.1111111111111111*nz*theta - phi*(0.1111111111111111 - 0.5*uz))/tauM;
            //........................................................................

            //Update velocity on device
    		Vel[0*Np+n] = ux;
    		Vel[1*Np+n] = uy;
    		Vel[2*Np+n] = uz;
            //Update pressure on device
            Pressure[n] = p;
            //Update chemical potential on device
            mu_phi[n] = chem;
            //Update color gradient on device
    		ColorGrad[0*Np+n] = nx;
    		ColorGrad[1*Np+n] = ny;
    		ColorGrad[2*Np+n] = nz;

    	}
	}
}

//__global__ void dvc_ScaLBL_D3Q19_AAodd_FreeLeeModel_Combined(int *neighborList, int *Map, double *dist, double *hq, double *Den,	double *Phi, double *mu_phi, double *Vel, double *Pressure, double *ColorGrad, 
//        double rhoA, double rhoB, double tauA, double tauB, double tauM, double kappa, double beta, double W, double Fx, double Fy, double Fz, 
//        int strideY, int strideZ, int start, int finish, int Np){
//
//	int n,nn,nn2x,ijk;
//	int nr1,nr2,nr3,nr4,nr5,nr6,nr7,nr8,nr9,nr10,nr11,nr12,nr13,nr14,nr15,nr16,nr17,nr18;
//    double ux,uy,uz;//fluid velocity 
//    double p;//pressure
//    double chem;//chemical potential
//    double phi; //phase field
//    double rho0;//fluid density
//	// distribution functions
//	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
//	double m0,m3,m5,m7;
//	double mm1,mm2,mm4,mm6,mm8,mm9,mm10,mm11,mm12,mm13,mm14,mm15,mm16,mm17,mm18;
//	double mm3,mm5,mm7;
//    double feq0,feq1,feq2,feq3,feq4,feq5,feq6,feq7,feq8,feq9,feq10,feq11,feq12,feq13,feq14,feq15,feq16,feq17,feq18;
//    double nx,ny,nz;//normal color gradient
//    double mgx,mgy,mgz;//mixed gradient reaching secondary neighbor
//
//    double h0,h1,h2,h3,h4,h5,h6;//distributions for LB phase field
//	double tau;//position dependent LB relaxation time for fluid
//    double C,theta;
//     double M = 2.0/9.0*(tauM-0.5);//diffusivity (or mobility) for the phase field D3Q7
//    double phi_temp;
//
//    int S = Np/NBLOCKS/NTHREADS + 1;
//    for (int s=0; s<S; s++){
//    	//........Get 1-D index for this thread....................
//    	n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
//
//    	if ( n<finish ){
//    		rho0 = Den[n];//load density
//
//    		// Get the 1D index based on regular data layout
//    		ijk = Map[n];
//            phi = Phi[ijk];// load phase field
//            phi_temp = phi;
//            if (phi>1.f) phi_temp=1.0;
//            if (phi<-1.f) phi_temp=-1.0;
//
//    		// local relaxation time
//    		tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);
//
//    		//					COMPUTE THE COLOR GRADIENT
//    		//........................................................................
//    		//.................Read Phase Indicator Values............................
//    		//........................................................................
//    		nn = ijk-1;							// neighbor index (get convention)
//    		m1 = Phi[nn];						// get neighbor for phi - 1
//    		//........................................................................
//    		nn = ijk+1;							// neighbor index (get convention)
//    		m2 = Phi[nn];						// get neighbor for phi - 2
//    		//........................................................................
//    		nn = ijk-strideY;							// neighbor index (get convention)
//    		m3 = Phi[nn];					// get neighbor for phi - 3
//    		//........................................................................
//    		nn = ijk+strideY;							// neighbor index (get convention)
//    		m4 = Phi[nn];					// get neighbor for phi - 4
//    		//........................................................................
//    		nn = ijk-strideZ;						// neighbor index (get convention)
//    		m5 = Phi[nn];					// get neighbor for phi - 5
//    		//........................................................................
//    		nn = ijk+strideZ;						// neighbor index (get convention)
//    		m6 = Phi[nn];					// get neighbor for phi - 6
//    		//........................................................................
//    		nn = ijk-strideY-1;						// neighbor index (get convention)
//    		m7 = Phi[nn];					// get neighbor for phi - 7
//    		//........................................................................
//    		nn = ijk+strideY+1;						// neighbor index (get convention)
//    		m8 = Phi[nn];					// get neighbor for phi - 8
//    		//........................................................................
//    		nn = ijk+strideY-1;						// neighbor index (get convention)
//    		m9 = Phi[nn];					// get neighbor for phi - 9
//    		//........................................................................
//    		nn = ijk-strideY+1;						// neighbor index (get convention)
//    		m10 = Phi[nn];					// get neighbor for phi - 10
//    		//........................................................................
//    		nn = ijk-strideZ-1;						// neighbor index (get convention)
//    		m11 = Phi[nn];					// get neighbor for phi - 11
//    		//........................................................................
//    		nn = ijk+strideZ+1;						// neighbor index (get convention)
//    		m12 = Phi[nn];					// get neighbor for phi - 12
//    		//........................................................................
//    		nn = ijk+strideZ-1;						// neighbor index (get convention)
//    		m13 = Phi[nn];					// get neighbor for phi - 13
//    		//........................................................................
//    		nn = ijk-strideZ+1;						// neighbor index (get convention)
//    		m14 = Phi[nn];					// get neighbor for phi - 14
//    		//........................................................................
//    		nn = ijk-strideZ-strideY;					// neighbor index (get convention)
//    		m15 = Phi[nn];					// get neighbor for phi - 15
//    		//........................................................................
//    		nn = ijk+strideZ+strideY;					// neighbor index (get convention)
//    		m16 = Phi[nn];					// get neighbor for phi - 16
//    		//........................................................................
//    		nn = ijk+strideZ-strideY;					// neighbor index (get convention)
//    		m17 = Phi[nn];					// get neighbor for phi - 17
//    		//........................................................................
//    		nn = ijk-strideZ+strideY;					// neighbor index (get convention)
//    		m18 = Phi[nn];					// get neighbor for phi - 18
//
//            // compute mixed difference (Eq.30, A.Fukhari et al. JCP 315(2016) 434-457)
//    		//........................................................................
//    		nn2x = ijk-2;							// neighbor index (get convention)
//    		mm1 = Phi[nn2x];						// get neighbor for phi - 1
//            mm1 = 0.25*(-mm1+5.0*m1-3.0*phi-m2);
//    		//........................................................................
//    		nn2x = ijk+2;							// neighbor index (get convention)
//    		mm2 = Phi[nn2x];						// get neighbor for phi - 2
//            mm2 = 0.25*(-mm2+5.0*m2-3.0*phi-m1);
//    		//........................................................................
//    		nn2x = ijk-strideY*2;							// neighbor index (get convention)
//    		mm3 = Phi[nn2x];					// get neighbor for phi - 3
//            mm3 = 0.25*(-mm3+5.0*m3-3.0*phi-m4);
//    		//........................................................................
//    		nn2x = ijk+strideY*2;							// neighbor index (get convention)
//    		mm4 = Phi[nn2x];					// get neighbor for phi - 4
//            mm4 = 0.25*(-mm4+5.0*m4-3.0*phi-m3);
//    		//........................................................................
//    		nn2x = ijk-strideZ*2;						// neighbor index (get convention)
//    		mm5 = Phi[nn2x];					// get neighbor for phi - 5
//            mm5 = 0.25*(-mm5+5.0*m5-3.0*phi-m6);
//    		//........................................................................
//    		nn2x = ijk+strideZ*2;						// neighbor index (get convention)
//    		mm6 = Phi[nn2x];					// get neighbor for phi - 6
//            mm6 = 0.25*(-mm6+5.0*m6-3.0*phi-m5);
//    		//........................................................................
//    		nn2x = ijk-strideY*2-2;						// neighbor index (get convention)
//    		mm7 = Phi[nn2x];					// get neighbor for phi - 7
//            mm7 = 0.25*(-mm7+5.0*m7-3.0*phi-m8);
//    		//........................................................................
//    		nn2x = ijk+strideY*2+2;						// neighbor index (get convention)
//    		mm8 = Phi[nn2x];					// get neighbor for phi - 8
//            mm8 = 0.25*(-mm8+5.0*m8-3.0*phi-m7);
//    		//........................................................................
//    		nn2x = ijk+strideY*2-2;						// neighbor index (get convention)
//    		mm9 = Phi[nn2x];					// get neighbor for phi - 9
//            mm9 = 0.25*(-mm9+5.0*m9-3.0*phi-m10);
//    		//........................................................................
//    		nn2x = ijk-strideY*2+2;						// neighbor index (get convention)
//    		mm10 = Phi[nn2x];					// get neighbor for phi - 10
//            mm10 = 0.25*(-mm10+5.0*m10-3.0*phi-m9);
//    		//........................................................................
//    		nn2x = ijk-strideZ*2-2;						// neighbor index (get convention)
//    		mm11 = Phi[nn2x];					// get neighbor for phi - 11
//            mm11 = 0.25*(-mm11+5.0*m11-3.0*phi-m12);
//    		//........................................................................
//    		nn2x = ijk+strideZ*2+2;						// neighbor index (get convention)
//    		mm12 = Phi[nn2x];					// get neighbor for phi - 12
//            mm12 = 0.25*(-mm12+5.0*m12-3.0*phi-m11);
//    		//........................................................................
//    		nn2x = ijk+strideZ*2-2;						// neighbor index (get convention)
//    		mm13 = Phi[nn2x];					// get neighbor for phi - 13
//            mm13 = 0.25*(-mm13+5.0*m13-3.0*phi-m14);
//    		//........................................................................
//    		nn2x = ijk-strideZ*2+2;						// neighbor index (get convention)
//    		mm14 = Phi[nn2x];					// get neighbor for phi - 14
//            mm14 = 0.25*(-mm14+5.0*m14-3.0*phi-m13);
//    		//........................................................................
//    		nn2x = ijk-strideZ*2-strideY*2;					// neighbor index (get convention)
//    		mm15 = Phi[nn2x];					// get neighbor for phi - 15
//            mm15 = 0.25*(-mm15+5.0*m15-3.0*phi-m16);
//    		//........................................................................
//    		nn2x = ijk+strideZ*2+strideY*2;					// neighbor index (get convention)
//    		mm16 = Phi[nn2x];					// get neighbor for phi - 16
//            mm16 = 0.25*(-mm16+5.0*m16-3.0*phi-m15);
//    		//........................................................................
//    		nn2x = ijk+strideZ*2-strideY*2;					// neighbor index (get convention)
//    		mm17 = Phi[nn2x];					// get neighbor for phi - 17
//            mm17 = 0.25*(-mm17+5.0*m17-3.0*phi-m18);
//    		//........................................................................
//    		nn2x = ijk-strideZ*2+strideY*2;					// neighbor index (get convention)
//    		mm18 = Phi[nn2x];					// get neighbor for phi - 18
//            mm18 = 0.25*(-mm18+5.0*m18-3.0*phi-m17);
//
//
//    		//............Compute the Color Gradient...................................
//    		nx = -3.0*1.0/18.0*(m1-m2+0.5*(m7-m8+m9-m10+m11-m12+m13-m14));
//    		ny = -3.0*1.0/18.0*(m3-m4+0.5*(m7-m8-m9+m10+m15-m16+m17-m18));
//    		nz = -3.0*1.0/18.0*(m5-m6+0.5*(m11-m12-m13+m14+m15-m16-m17+m18));
//    		//............Compute the Chemical Potential...............................
//            //chem = 2.0*3.0/18.0*(m1+m2+m3+m4+m5+m6-6*phi+0.5*(m7+m8+m9+m10+m11+m12+m13+m14+m15+m16+m17+m18-12*phi));//intermediate var, i.e. the laplacian
//            //chem = 4.0*beta*phi*(phi+1.0)*(phi-1.0)-kappa*chem;
//            chem = 2.0*3.0/18.0*(m1+m2+m3+m4+m5+m6-6*phi_temp+0.5*(m7+m8+m9+m10+m11+m12+m13+m14+m15+m16+m17+m18-12*phi_temp));//intermediate var, i.e. the laplacian
//            chem = 4.0*beta*phi_temp*(phi_temp+1.0)*(phi_temp-1.0)-kappa*chem;
//    		//............Compute the Mixed Gradient...................................
//    		mgx = -3.0*1.0/18.0*(mm1-mm2+0.5*(mm7-mm8+mm9-mm10+mm11-mm12+mm13-mm14));
//    		mgy = -3.0*1.0/18.0*(mm3-mm4+0.5*(mm7-mm8-mm9+mm10+mm15-mm16+mm17-mm18));
//    		mgz = -3.0*1.0/18.0*(mm5-mm6+0.5*(mm11-mm12-mm13+mm14+mm15-mm16-mm17+mm18));
//    		
//            //de-noise color gradient and mixed gradient
//            C = sqrt(nx*nx+ny*ny+nz*nz);
//            if (C<1.0e-12) nx=ny=nz=0.0;
//            double mg_mag = sqrt(mgx*mgx+mgy*mgy+mgz*mgz);
//            if (mg_mag<1.0e-12) mgx=mgy=mgz=0.0;
//            //maybe you can also de-noise chemical potential ? within the bulk phase chem should be zero
//            if (fabs(chem)<1.0e-12) chem=0.0;
//
//    		// q=0
//    		m0 = dist[n];
//    		// q=1
//    		nr1 = neighborList[n]; // neighbor 2 ( > 10Np => odd part of dist)
//    		m1 = dist[nr1]; // reading the f1 data into register fq
//
//    		nr2 = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
//    		m2 = dist[nr2];  // reading the f2 data into register fq
//
//    		// q=3
//    		nr3 = neighborList[n+2*Np]; // neighbor 4
//    		m3 = dist[nr3];
//
//    		// q = 4
//    		nr4 = neighborList[n+3*Np]; // neighbor 3
//    		m4 = dist[nr4];
//
//    		// q=5
//    		nr5 = neighborList[n+4*Np];
//    		m5 = dist[nr5];
//
//    		// q = 6
//    		nr6 = neighborList[n+5*Np];
//    		m6 = dist[nr6];
//    		
//    		// q=7
//    		nr7 = neighborList[n+6*Np];
//    		m7 = dist[nr7];
//
//    		// q = 8
//    		nr8 = neighborList[n+7*Np];
//    		m8 = dist[nr8];
//
//    		// q=9
//    		nr9 = neighborList[n+8*Np];
//    		m9 = dist[nr9];
//
//    		// q = 10
//    		nr10 = neighborList[n+9*Np];
//    		m10 = dist[nr10];
//
//    		// q=11
//    		nr11 = neighborList[n+10*Np];
//    		m11 = dist[nr11];
//
//    		// q=12
//    		nr12 = neighborList[n+11*Np];
//    		m12 = dist[nr12];
//
//    		// q=13
//    		nr13 = neighborList[n+12*Np];
//    		m13 = dist[nr13];
//
//    		// q=14
//    		nr14 = neighborList[n+13*Np];
//    		m14 = dist[nr14];
//
//    		// q=15
//    		nr15 = neighborList[n+14*Np];
//    		m15 = dist[nr15];
//
//    		// q=16
//    		nr16 = neighborList[n+15*Np];
//    		m16 = dist[nr16];
//
//    		// q=17
//    		nr17 = neighborList[n+16*Np];
//    		m17 = dist[nr17];
//
//    		// q=18
//    		nr18 = neighborList[n+17*Np];
//    		m18 = dist[nr18];
//
//            //compute fluid velocity
//            ux = 3.0/rho0*(m1-m2+m7-m8+m9-m10+m11-m12+m13-m14+0.5*(chem*nx+Fx)/3.0);
//            uy = 3.0/rho0*(m3-m4+m7-m8-m9+m10+m15-m16+m17-m18+0.5*(chem*ny+Fy)/3.0);
//            uz = 3.0/rho0*(m5-m6+m11-m12-m13+m14+m15-m16-m17+m18+0.5*(chem*nz+Fz)/3.0);
//            //compute pressure
//            p = (m0+m2+m1+m4+m3+m6+m5+m8+m7+m10+m9+m12+m11+m14+m13+m16+m15+m18+m17)
//                      +0.5*(rhoA-rhoB)/2.0/3.0*(ux*nx+uy*ny+uz*nz);
//
//            //compute equilibrium distributions
//                feq0 = 0.3333333333333333*p - 0.25*(Fx*ux + Fy*uy + Fz*uz)*(-0.6666666666666666 + ux*ux + uy*uy + uz*uz) - 
//         0.16666666666666666*rho0*(ux*ux + uy*uy + uz*uz) - 0.5*(-(nx*ux) - ny*uy - nz*uz)*
//          (-0.08333333333333333*(rhoA - rhoB)*(ux*ux + uy*uy + uz*uz) + chem*(0.3333333333333333 - 0.5*(ux*ux + uy*uy + uz*uz)));
//                feq1 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-ux*ux + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*uz)) - 
//         0.125*(Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 
//           0.3333333333333333*(-2.*ux + ux*ux + uy*uy + uz*uz)) - 0.0625*(nx - nx*ux - ny*uy - nz*uz)*
//          (2*chem*ux*ux - 0.3333333333333333*((-rhoA + rhoB)*ux*ux + 2*chem*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
//           0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*uz)));
//                feq2 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-ux*ux + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*uz)) - 
//         0.125*(Fx + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 
//           0.3333333333333333*(2.*ux + ux*ux + uy*uy + uz*uz)) - 0.0625*(nx + nx*ux + ny*uy + nz*uz)*
//          (-2.*chem*ux*ux + 0.1111111111111111*(-4.*chem + rhoB*(-2.*ux - 1.*ux*ux - 1.*uy*uy - 1.*uz*uz) + 
//             rhoA*(2.*ux + ux*ux + uy*uy + uz*uz)) + 0.3333333333333333*((-1.*rhoA + rhoB)*ux*ux + 
//             chem*(4.*ux + 2.*ux*ux + 2.*uy*uy + 2.*uz*uz)));
//                feq3 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uy*uy + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*uz)) - 
//         0.125*(Fx*ux + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy + 
//           0.3333333333333333*(ux*ux - 2.*uy + uy*uy + uz*uz)) - 0.0625*(ny - nx*ux - ny*uy - nz*uz)*
//          (2*chem*uy*uy - 0.3333333333333333*((-rhoA + rhoB)*uy*uy + 2*chem*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
//           0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*uz)));
//                feq4 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uy*uy + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*uz)) - 
//         0.125*(Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy + 
//           0.3333333333333333*(ux*ux + 2.*uy + uy*uy + uz*uz)) - 0.0625*(ny + nx*ux + ny*uy + nz*uz)*
//          (-2.*chem*uy*uy + 0.1111111111111111*(-4.*chem + rhoB*(-1.*ux*ux - 2.*uy - 1.*uy*uy - 1.*uz*uz) + 
//             rhoA*(ux*ux + 2.*uy + uy*uy + uz*uz)) + 0.3333333333333333*((-1.*rhoA + rhoB)*uy*uy + 
//             chem*(2.*ux*ux + 4.*uy + 2.*uy*uy + 2.*uz*uz))); 
//                feq5 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uz*uz + 0.3333333333333333*(ux*ux + uy*uy + (-2 + uz)*uz)) - 
//         0.125*(Fx*ux + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*uz*uz + 
//           0.3333333333333333*(ux*ux + uy*uy + (-2. + uz)*uz)) - 0.0625*(nx*ux + ny*uy + nz*(-1. + uz))*
//          (-2.*chem*uz*uz + 0.1111111111111111*(-4.*chem + rhoB*(-1.*ux*ux - 1.*uy*uy + (2. - 1.*uz)*uz) + 
//             rhoA*(ux*ux + uy*uy + (-2. + uz)*uz)) + 0.3333333333333333*((-1.*rhoA + rhoB)*uz*uz + 
//             chem*(2.*ux*ux + 2.*uy*uy + uz*(-4. + 2.*uz)))); 
//                feq6 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uz*uz + 0.3333333333333333*(ux*ux + uy*uy + uz*(2 + uz))) - 
//         0.125*(Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*uz*uz + 
//           0.3333333333333333*(ux*ux + uy*uy + uz*(2. + uz))) - 0.0625*(nz + nx*ux + ny*uy + nz*uz)*
//          (-2.*chem*uz*uz + 0.1111111111111111*(-4.*chem + rhoB*(-1.*ux*ux - 1.*uy*uy + (-2. - 1.*uz)*uz) + 
//             rhoA*(ux*ux + uy*uy + uz*(2. + uz))) + 0.3333333333333333*((-1.*rhoA + rhoB)*uz*uz + 
//             chem*(2.*ux*ux + 2.*uy*uy + uz*(4. + 2.*uz)))); 
//                feq7 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(ux + uy)*(ux + uy) + 0.3333333333333333*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) - 
//         0.0625*(Fx*(-1. + ux) + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uy - 1.*uy*uy + 
//           0.3333333333333333*(-2.*ux + ux*ux - 2.*uy + uy*uy + uz*uz)) - 0.03125*(nx + ny - nx*ux - ny*uy - nz*uz)*
//          (2*chem*(ux + uy)*(ux + uy) + 0.3333333333333333*((rhoA - rhoB)*(ux + uy)*(ux + uy) - 2*chem*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
//           0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz))); 
//                feq8 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(ux + uy)*(ux + uy) + 0.3333333333333333*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) - 
//         0.0625*(Fx + Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uy - 1.*uy*uy + 
//           0.3333333333333333*(2.*ux + ux*ux + 2.*uy + uy*uy + uz*uz)) - 0.03125*(-(nx*(1 + ux)) - ny*(1 + uy) - nz*uz)*
//          (2*chem*(ux + uy)*(ux + uy) - 0.3333333333333333*(-((rhoA - rhoB)*(ux + uy)*(ux + uy)) + 
//             2*chem*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 0.1111111111111111*
//            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz))); 
//                feq9 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(ux - uy)*(ux - uy) + 0.3333333333333333*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) - 
//         0.0625*(Fy + Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uy - 1.*uy*uy + 
//           0.3333333333333333*(-2.*ux + ux*ux + 2.*uy + uy*uy + uz*uz)) - 0.03125*(nx - nx*ux - ny*(1 + uy) - nz*uz)*
//          (2*chem*(ux - uy)*(ux - uy) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uy)*(ux - uy)) + 
//             2*chem*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 0.1111111111111111*
//            (4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz))); 
//                feq10 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(ux - uy)*(ux - uy) + 0.3333333333333333*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) - 
//         0.0625*(Fx*(1 + ux) + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uy - 1.*uy*uy + 
//           0.3333333333333333*(2.*ux + ux*ux - 2.*uy + uy*uy + uz*uz)) - 0.03125*(ny - nx*(1 + ux) - ny*uy - nz*uz)*
//          (2*chem*(ux - uy)*(ux - uy) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uy)*(ux - uy)) + 
//             2*chem*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 0.1111111111111111*
//            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz))); 
//                feq11 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(ux + uz)*(ux + uz) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) - 
//         0.0625*(Fx*(-1. + ux) + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uz - 1.*uz*uz + 
//           0.3333333333333333*(-2.*ux + ux*ux + uy*uy + (-2. + uz)*uz)) - 0.03125*(nx + nz - nx*ux - ny*uy - nz*uz)*
//          (2*chem*(ux + uz)*(ux + uz) + 0.3333333333333333*((rhoA - rhoB)*(ux + uz)*(ux + uz) - 2*chem*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
//           0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz))); 
//                feq12 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(ux + uz)*(ux + uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) - 
//         0.0625*(Fx + Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uz - 1.*uz*uz + 
//           0.3333333333333333*(2.*ux + ux*ux + uy*uy + uz*(2. + uz))) - 0.03125*(-(nx*(1 + ux)) - ny*uy - nz*(1 + uz))*
//          (2*chem*(ux + uz)*(ux + uz) - 0.3333333333333333*(-((rhoA - rhoB)*(ux + uz)*(ux + uz)) + 
//             2*chem*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
//            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))); 
//                feq13 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(ux - uz)*(ux - uz) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) - 
//         0.0625*(Fz + Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uz - 1.*uz*uz + 
//           0.3333333333333333*(-2.*ux + ux*ux + uy*uy + uz*(2. + uz))) - 0.03125*(nx - nx*ux - ny*uy - nz*(1 + uz))*
//          (2*chem*(ux - uz)*(ux - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uz)*(ux - uz)) + 
//             2*chem*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
//            (4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*(2 + uz)))); 
//                feq14 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(ux - uz)*(ux - uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) - 
//         0.0625*(Fx*(1 + ux) + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uz - 1.*uz*uz + 
//           0.3333333333333333*(2.*ux + ux*ux + uy*uy + (-2. + uz)*uz)) - 0.03125*(nz - nx*(1 + ux) - ny*uy - nz*uz)*
//          (2*chem*(ux - uz)*(ux - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uz)*(ux - uz)) + 
//             2*chem*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 0.1111111111111111*
//            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz))); 
//                feq15 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) - 
//         0.0625*(Fx*ux + Fy*(-1. + uy) + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*uy*uy - 2.*uy*uz - 1.*uz*uz + 
//           0.3333333333333333*(ux*ux - 2.*uy + uy*uy + (-2. + uz)*uz)) - 0.03125*(ny + nz - nx*ux - ny*uy - nz*uz)*
//          (2*chem*(uy + uz)*(uy + uz) + 0.3333333333333333*((rhoA - rhoB)*(uy + uz)*(uy + uz) - 2*chem*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) + 
//           0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))); 
//                feq16 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) - 
//         0.0625*(Fy + Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy - 2.*uy*uz - 1.*uz*uz + 
//           0.3333333333333333*(ux*ux + 2.*uy + uy*uy + uz*(2. + uz))) - 0.03125*(-(nx*ux) - ny*(1 + uy) - nz*(1 + uz))*
//          (2*chem*(uy + uz)*(uy + uz) - 0.3333333333333333*(-((rhoA - rhoB)*(uy + uz)*(uy + uz)) + 
//             2*chem*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
//            (4*chem - (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))); 
//                feq17 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) - 
//         0.0625*(Fz + Fx*ux + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy + 2.*uy*uz - 1.*uz*uz + 
//           0.3333333333333333*(ux*ux - 2.*uy + uy*uy + uz*(2. + uz))) - 0.03125*(ny - nx*ux - ny*uy - nz*(1 + uz))*
//          (2*chem*(uy - uz)*(uy - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(uy - uz)*(uy - uz)) + 
//             2*chem*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
//            (4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))); 
//                feq18 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) - 
//         0.0625*(Fx*ux + Fy*(1 + uy) + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*uy*uy + 2.*uy*uz - 1.*uz*uz + 
//           0.3333333333333333*(ux*ux + 2.*uy + uy*uy + (-2. + uz)*uz)) - 0.03125*(nz - nx*ux - ny*(1 + uy) - nz*uz)*
//          (2*chem*(uy - uz)*(uy - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(uy - uz)*(uy - uz)) + 
//             2*chem*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 0.1111111111111111*
//            (4*chem - (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz))); 
//
//            //------------------------------------------------- BCK collison ------------------------------------------------------------//
//    		// q=0
//    		dist[n] = m0 - (m0-feq0)/tau + 0.25*(2*(Fx*ux + Fy*uy + Fz*uz)*(-0.6666666666666666 + ux*ux + uy*uy + uz*uz) + 
//         (mgx*ux + mgy*uy + mgz*uz)*(2*chem*(ux*ux + uy*uy + uz*uz) + 
//            0.3333333333333333*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + uz*uz)))); 
//
//    		// q = 1
//    		dist[nr2] = m1 - (m1-feq1)/tau + 0.125*(2*(Fx*(-1 + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - ux*ux + 
//            0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
//         (mgx*(-1 + ux) + mgy*uy + mgz*uz)*(-2*chem*(ux*ux) + 
//            0.3333333333333333*((-rhoA + rhoB)*(ux*ux) + 2*chem*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*uz))));
//
//    		// q=2
//    		dist[nr1] = m2 - (m2-feq2)/tau + 0.125*(2*(Fx + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - ux*ux + 
//            0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*uz)) + 
//         (mgx + mgx*ux + mgy*uy + mgz*uz)*(-2*chem*(ux*ux) + 
//            0.3333333333333333*((-rhoA + rhoB)*(ux*ux) + 2*chem*(2*ux + ux*ux + uy*uy + uz*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*uz))));
//
//    		// q = 3
//    		dist[nr4] = m3 - (m3-feq3)/tau + 0.125*(2*(Fx*ux + Fy*(-1 + uy) + Fz*uz)*(-0.2222222222222222 - uy*uy + 
//            0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
//         (mgx*ux + mgy*(-1 + uy) + mgz*uz)*(-2*chem*(uy*uy) + 
//            0.3333333333333333*((-rhoA + rhoB)*(uy*uy) + 2*chem*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*uz))));
//
//    		// q = 4
//    		dist[nr3] = m4 - (m4-feq4)/tau + 0.125*(2*(Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - uy*uy + 
//            0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*uz)) + 
//         (mgy + mgx*ux + mgy*uy + mgz*uz)*(-2*chem*(uy*uy) + 
//            0.3333333333333333*((-rhoA + rhoB)*(uy*uy) + 2*chem*(ux*ux + 2*uy + uy*uy + uz*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*uz))));
//
//    		// q = 5
//    		dist[nr6] = m5 - (m5-feq5)/tau + 0.125*(2*(Fx*ux + Fy*uy + Fz*(-1 + uz))*(-0.2222222222222222 - uz*uz + 
//            0.3333333333333333*(ux*ux + uy*uy + (-2 + uz)*uz)) + 
//         (mgx*ux + mgy*uy + mgz*(-1 + uz))*(-2*chem*(uz*uz) + 
//            0.3333333333333333*((-rhoA + rhoB)*(uz*uz) + 2*chem*(ux*ux + uy*uy + (-2 + uz)*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + (-2 + uz)*uz))));
//
//    		// q = 6
//    		dist[nr5] = m6 - (m6-feq6)/tau + 0.125*(2*(Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - uz*uz + 
//            0.3333333333333333*(ux*ux + uy*uy + uz*(2 + uz))) + 
//         (mgz + mgx*ux + mgy*uy + mgz*uz)*(-2*chem*(uz*uz) + 
//            0.3333333333333333*((-rhoA + rhoB)*(uz*uz) + 2*chem*(ux*ux + uy*uy + uz*(2 + uz))) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + uz*(2 + uz)))));
//
//    		// q = 7
//    		dist[nr8] = m7 - (m7-feq7)/tau + 0.0625*(-2*(Fx*(-1 + ux) + Fy*(-1 + uy) + Fz*uz)*
//          (0.2222222222222222 + (ux + uy)*(ux + uy) - 
//            0.3333333333333333*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
//         (mgx*(-1 + ux) + mgy*(-1 + uy) + mgz*uz)*
//          (-2*chem*((ux + uy)*(ux + uy)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((ux + uy)*(ux + uy))) + 2*chem*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz))));
//
//    		// q = 8
//    		dist[nr7] = m8 - (m8-feq8)/tau + 0.0625*(2*(Fx + Fy + Fx*ux + Fy*uy + Fz*uz)*
//          (-0.2222222222222222 - (ux + uy)*(ux + uy) + 
//            0.3333333333333333*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
//         (mgx + mgy + mgx*ux + mgy*uy + mgz*uz)*
//          (-2*chem*((ux + uy)*(ux + uy)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((ux + uy)*(ux + uy))) + 2*chem*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz))));
//
//    		// q = 9
//    		dist[nr10] = m9 - (m9-feq9)/tau + 0.0625*(2*(Fy + Fx*(-1 + ux) + Fy*uy + Fz*uz)*
//          (-0.2222222222222222 - (ux - uy)*(ux - uy) + 
//            0.3333333333333333*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
//         (mgy + mgx*(-1 + ux) + mgy*uy + mgz*uz)*
//          (-2*chem*((ux - uy)*(ux - uy)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((ux - uy)*(ux - uy))) + 2*chem*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz))));
//
//    		// q = 10
//    		dist[nr9] = m10 - (m10-feq10)/tau + 0.0625*(2*(Fx*(1 + ux) + Fy*(-1 + uy) + Fz*uz)*
//          (-0.2222222222222222 - (ux - uy)*(ux - uy) + 
//            0.3333333333333333*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
//         (mgx*(1 + ux) + mgy*(-1 + uy) + mgz*uz)*
//          (-2*chem*((ux - uy)*(ux - uy)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((ux - uy)*(ux - uy))) + 2*chem*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz))));
//
//    		// q = 11
//    		dist[nr12] = m11 - (m11-feq11)/tau + 0.0625*(-2*(Fx*(-1 + ux) + Fy*uy + Fz*(-1 + uz))*
//          (0.2222222222222222 + (ux + uz)*(ux + uz) - 
//            0.3333333333333333*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
//         (mgx*(-1 + ux) + mgy*uy + mgz*(-1 + uz))*
//          (-2*chem*((ux + uz)*(ux + uz)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((ux + uz)*(ux + uz))) + 2*chem*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz))));
//
//    		// q = 12
//    		dist[nr11] = m12 - (m12-feq12)/tau + 0.0625*(2*(Fx + Fz + Fx*ux + Fy*uy + Fz*uz)*
//          (-0.2222222222222222 - (ux + uz)*(ux + uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))
//           + (mgx + mgz + mgx*ux + mgy*uy + mgz*uz)*
//          (-2*chem*((ux + uz)*(ux + uz)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((ux + uz)*(ux + uz))) + 2*chem*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))));
//
//    		// q = 13
//    		dist[nr14] = m13 - (m13-feq13)/tau + 0.0625*(2*(Fz + Fx*(-1 + ux) + Fy*uy + Fz*uz)*
//          (-0.2222222222222222 - (ux - uz)*(ux - uz) + 
//            0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
//         (mgz + mgx*(-1 + ux) + mgy*uy + mgz*uz)*
//          (-2*chem*((ux - uz)*(ux - uz)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((ux - uz)*(ux - uz))) + 2*chem*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*(2 + uz)))));
//
//    		// q= 14
//    		dist[nr13] = m14 - (m14-feq14)/tau + 0.0625*(2*(Fx*(1 + ux) + Fy*uy + Fz*(-1 + uz))*
//          (-0.2222222222222222 - (ux - uz)*(ux - uz) + 
//            0.3333333333333333*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
//         (mgx*(1 + ux) + mgy*uy + mgz*(-1 + uz))*
//          (-2*chem*((ux - uz)*(ux - uz)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((ux - uz)*(ux - uz))) + 2*chem*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz))));
//
//    		// q = 15
//    		dist[nr16] = m15 - (m15-feq15)/tau + 0.0625*(-2*(Fx*ux + Fy*(-1 + uy) + Fz*(-1 + uz))*
//          (0.2222222222222222 + (uy + uz)*(uy + uz) - 0.3333333333333333*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))
//           + (mgx*ux + mgy*(-1 + uy) + mgz*(-1 + uz))*
//          (-2*chem*((uy + uz)*(uy + uz)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((uy + uz)*(uy + uz))) + 2*chem*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))));
//
//    		// q = 16
//    		dist[nr15] = m16 - (m16-feq16)/tau + 0.0625*(2*(Fy + Fz + Fx*ux + Fy*uy + Fz*uz)*
//          (-0.2222222222222222 - (uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))
//           + (mgy + mgz + mgx*ux + mgy*uy + mgz*uz)*
//          (-2*chem*((uy + uz)*(uy + uz)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((uy + uz)*(uy + uz))) + 2*chem*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))));
//
//    		// q = 17
//    		dist[nr18] = m17 - (m17-feq17)/tau + 0.0625*(2*(Fz + Fx*ux + Fy*(-1 + uy) + Fz*uz)*
//          (-0.2222222222222222 - (uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))
//           + (mgz + mgx*ux + mgy*(-1 + uy) + mgz*uz)*
//          (-2*chem*((uy - uz)*(uy - uz)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((uy - uz)*(uy - uz))) + 2*chem*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))));
//
//    		// q = 18
//    		dist[nr17] = m18 - (m18-feq18)/tau + 0.0625*(2*(Fx*ux + Fy*(1 + uy) + Fz*(-1 + uz))*
//          (-0.2222222222222222 - (uy - uz)*(uy - uz) + 
//            0.3333333333333333*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 
//         (mgx*ux + mgy*(1 + uy) + mgz*(-1 + uz))*
//          (-2*chem*((uy - uz)*(uy - uz)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((uy - uz)*(uy - uz))) + 2*chem*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz))));
//            //----------------------------------------------------------------------------------------------------------------------------------------//
//
//            // ----------------------------- compute phase field evolution ----------------------------------------
//            //Normalize the Color Gradient
//            C = sqrt(nx*nx+ny*ny+nz*nz);
//            double ColorMag = C;
//            if (C==0.0) ColorMag=1.0;
//            nx = nx/ColorMag;
//            ny = ny/ColorMag;
//            nz = nz/ColorMag;		
//            //compute surface tension-related parameter
//            //theta = 4.5*M*2.0*(1-phi*phi)/W;
//            theta = 4.5*M*2.0*(1-phi_temp*phi_temp)/W;
//
//            //load distributions of phase field
//            //q=0
//            h0 = hq[n];
//            //q=1
//            h1 = hq[nr1]; 
//
//            //q=2
//            h2 = hq[nr2];  
//
//            //q=3
//            h3 = hq[nr3];
//
//            //q=4
//            h4 = hq[nr4];
//
//            //q=5
//            h5 = hq[nr5];
//
//            //q=6
//            h6 = hq[nr6];
//
//            //-------------------------------- BGK collison for phase field ---------------------------------//
//            // q = 0
//            hq[n] = h0 - (h0 - 0.3333333333333333*phi)/tauM;
//
//            // q = 1
//            hq[nr2] = h1 - (h1 - 0.1111111111111111*nx*theta - phi*(0.1111111111111111 + 0.5*ux))/tauM;
//
//            // q = 2
//            hq[nr1] = h2 - (h2 + 0.1111111111111111*nx*theta - phi*(0.1111111111111111 - 0.5*ux))/tauM;
//
//            // q = 3
//            hq[nr4] = h3 - (h3 - 0.1111111111111111*ny*theta - phi*(0.1111111111111111 + 0.5*uy))/tauM;
//
//            // q = 4
//            hq[nr3] = h4 - (h4 + 0.1111111111111111*ny*theta - phi*(0.1111111111111111 - 0.5*uy))/tauM;
//
//            // q = 5
//            hq[nr6] = h5 - (h5 - 0.1111111111111111*nz*theta - phi*(0.1111111111111111 + 0.5*uz))/tauM;
//
//            // q = 6
//            hq[nr5] = h6 - (h6 + 0.1111111111111111*nz*theta - phi*(0.1111111111111111 - 0.5*uz))/tauM;
//            //........................................................................
//
//            //Update velocity on device
//    		Vel[0*Np+n] = ux;
//    		Vel[1*Np+n] = uy;
//    		Vel[2*Np+n] = uz;
//            //Update pressure on device
//            Pressure[n] = p;
//            //Update chemical potential on device
//            mu_phi[n] = chem;
//            //Update color gradient on device
//    		ColorGrad[0*Np+n] = nx;
//    		ColorGrad[1*Np+n] = ny;
//    		ColorGrad[2*Np+n] = nz;
//    	}
//    }
//}

//__global__ void dvc_ScaLBL_D3Q19_AAeven_FreeLeeModel_Combined(int *Map, double *dist, double *hq, double *Den,	double *Phi, double *mu_phi, double *Vel, double *Pressure, double *ColorGrad,
//        double rhoA, double rhoB, double tauA, double tauB, double tauM, double kappa, double beta, double W, double Fx, double Fy, double Fz, 
//        int strideY, int strideZ, int start, int finish, int Np){
//
//	int n,nn,nn2x,ijk;
//    double ux,uy,uz;//fluid velocity 
//    double p;//pressure
//    double chem;//chemical potential
//    double phi; //phase field
//    double rho0;//fluid density
//	// distribution functions
//	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
//	double m0,m3,m5,m7;
//	double mm1,mm2,mm4,mm6,mm8,mm9,mm10,mm11,mm12,mm13,mm14,mm15,mm16,mm17,mm18;
//	double mm3,mm5,mm7;
//    double feq0,feq1,feq2,feq3,feq4,feq5,feq6,feq7,feq8,feq9,feq10,feq11,feq12,feq13,feq14,feq15,feq16,feq17,feq18;
//    double nx,ny,nz;//normal color gradient
//    double mgx,mgy,mgz;//mixed gradient reaching secondary neighbor
//
//    double h0,h1,h2,h3,h4,h5,h6;//distributions for LB phase field
//	double tau;//position dependent LB relaxation time for fluid
//    double C,theta;
//    double M = 2.0/9.0*(tauM-0.5);//diffusivity (or mobility) for the phase field D3Q7
//    double phi_temp;
//
//    int S = Np/NBLOCKS/NTHREADS + 1;
//    for (int s=0; s<S; s++){
//    	//........Get 1-D index for this thread....................
//    	n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
//
//    	if ( n<finish ){
//    		rho0 = Den[n];//load density
//
//    		// Get the 1D index based on regular data layout
//    		ijk = Map[n];
//            phi = Phi[ijk];// load phase field
//            phi_temp = phi;
//            if (phi>1.f) phi_temp=1.0;
//            if (phi<-1.f) phi_temp=-1.0;
//
//    		// local relaxation time
//    		tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);
//
//    		//					COMPUTE THE COLOR GRADIENT
//    		//........................................................................
//    		//.................Read Phase Indicator Values............................
//    		//........................................................................
//    		nn = ijk-1;							// neighbor index (get convention)
//    		m1 = Phi[nn];						// get neighbor for phi - 1
//    		//........................................................................
//    		nn = ijk+1;							// neighbor index (get convention)
//    		m2 = Phi[nn];						// get neighbor for phi - 2
//    		//........................................................................
//    		nn = ijk-strideY;							// neighbor index (get convention)
//    		m3 = Phi[nn];					// get neighbor for phi - 3
//    		//........................................................................
//    		nn = ijk+strideY;							// neighbor index (get convention)
//    		m4 = Phi[nn];					// get neighbor for phi - 4
//    		//........................................................................
//    		nn = ijk-strideZ;						// neighbor index (get convention)
//    		m5 = Phi[nn];					// get neighbor for phi - 5
//    		//........................................................................
//    		nn = ijk+strideZ;						// neighbor index (get convention)
//    		m6 = Phi[nn];					// get neighbor for phi - 6
//    		//........................................................................
//    		nn = ijk-strideY-1;						// neighbor index (get convention)
//    		m7 = Phi[nn];					// get neighbor for phi - 7
//    		//........................................................................
//    		nn = ijk+strideY+1;						// neighbor index (get convention)
//    		m8 = Phi[nn];					// get neighbor for phi - 8
//    		//........................................................................
//    		nn = ijk+strideY-1;						// neighbor index (get convention)
//    		m9 = Phi[nn];					// get neighbor for phi - 9
//    		//........................................................................
//    		nn = ijk-strideY+1;						// neighbor index (get convention)
//    		m10 = Phi[nn];					// get neighbor for phi - 10
//    		//........................................................................
//    		nn = ijk-strideZ-1;						// neighbor index (get convention)
//    		m11 = Phi[nn];					// get neighbor for phi - 11
//    		//........................................................................
//    		nn = ijk+strideZ+1;						// neighbor index (get convention)
//    		m12 = Phi[nn];					// get neighbor for phi - 12
//    		//........................................................................
//    		nn = ijk+strideZ-1;						// neighbor index (get convention)
//    		m13 = Phi[nn];					// get neighbor for phi - 13
//    		//........................................................................
//    		nn = ijk-strideZ+1;						// neighbor index (get convention)
//    		m14 = Phi[nn];					// get neighbor for phi - 14
//    		//........................................................................
//    		nn = ijk-strideZ-strideY;					// neighbor index (get convention)
//    		m15 = Phi[nn];					// get neighbor for phi - 15
//    		//........................................................................
//    		nn = ijk+strideZ+strideY;					// neighbor index (get convention)
//    		m16 = Phi[nn];					// get neighbor for phi - 16
//    		//........................................................................
//    		nn = ijk+strideZ-strideY;					// neighbor index (get convention)
//    		m17 = Phi[nn];					// get neighbor for phi - 17
//    		//........................................................................
//    		nn = ijk-strideZ+strideY;					// neighbor index (get convention)
//    		m18 = Phi[nn];					// get neighbor for phi - 18
//
//            // compute mixed difference (Eq.30, A.Fukhari et al. JCP 315(2016) 434-457)
//    		//........................................................................
//    		nn2x = ijk-2;							// neighbor index (get convention)
//    		mm1 = Phi[nn2x];						// get neighbor for phi - 1
//            mm1 = 0.25*(-mm1+5.0*m1-3.0*phi-m2);
//    		//........................................................................
//    		nn2x = ijk+2;							// neighbor index (get convention)
//    		mm2 = Phi[nn2x];						// get neighbor for phi - 2
//            mm2 = 0.25*(-mm2+5.0*m2-3.0*phi-m1);
//    		//........................................................................
//    		nn2x = ijk-strideY*2;							// neighbor index (get convention)
//    		mm3 = Phi[nn2x];					// get neighbor for phi - 3
//            mm3 = 0.25*(-mm3+5.0*m3-3.0*phi-m4);
//    		//........................................................................
//    		nn2x = ijk+strideY*2;							// neighbor index (get convention)
//    		mm4 = Phi[nn2x];					// get neighbor for phi - 4
//            mm4 = 0.25*(-mm4+5.0*m4-3.0*phi-m3);
//    		//........................................................................
//    		nn2x = ijk-strideZ*2;						// neighbor index (get convention)
//    		mm5 = Phi[nn2x];					// get neighbor for phi - 5
//            mm5 = 0.25*(-mm5+5.0*m5-3.0*phi-m6);
//    		//........................................................................
//    		nn2x = ijk+strideZ*2;						// neighbor index (get convention)
//    		mm6 = Phi[nn2x];					// get neighbor for phi - 6
//            mm6 = 0.25*(-mm6+5.0*m6-3.0*phi-m5);
//    		//........................................................................
//    		nn2x = ijk-strideY*2-2;						// neighbor index (get convention)
//    		mm7 = Phi[nn2x];					// get neighbor for phi - 7
//            mm7 = 0.25*(-mm7+5.0*m7-3.0*phi-m8);
//    		//........................................................................
//    		nn2x = ijk+strideY*2+2;						// neighbor index (get convention)
//    		mm8 = Phi[nn2x];					// get neighbor for phi - 8
//            mm8 = 0.25*(-mm8+5.0*m8-3.0*phi-m7);
//    		//........................................................................
//    		nn2x = ijk+strideY*2-2;						// neighbor index (get convention)
//    		mm9 = Phi[nn2x];					// get neighbor for phi - 9
//            mm9 = 0.25*(-mm9+5.0*m9-3.0*phi-m10);
//    		//........................................................................
//    		nn2x = ijk-strideY*2+2;						// neighbor index (get convention)
//    		mm10 = Phi[nn2x];					// get neighbor for phi - 10
//            mm10 = 0.25*(-mm10+5.0*m10-3.0*phi-m9);
//    		//........................................................................
//    		nn2x = ijk-strideZ*2-2;						// neighbor index (get convention)
//    		mm11 = Phi[nn2x];					// get neighbor for phi - 11
//            mm11 = 0.25*(-mm11+5.0*m11-3.0*phi-m12);
//    		//........................................................................
//    		nn2x = ijk+strideZ*2+2;						// neighbor index (get convention)
//    		mm12 = Phi[nn2x];					// get neighbor for phi - 12
//            mm12 = 0.25*(-mm12+5.0*m12-3.0*phi-m11);
//    		//........................................................................
//    		nn2x = ijk+strideZ*2-2;						// neighbor index (get convention)
//    		mm13 = Phi[nn2x];					// get neighbor for phi - 13
//            mm13 = 0.25*(-mm13+5.0*m13-3.0*phi-m14);
//    		//........................................................................
//    		nn2x = ijk-strideZ*2+2;						// neighbor index (get convention)
//    		mm14 = Phi[nn2x];					// get neighbor for phi - 14
//            mm14 = 0.25*(-mm14+5.0*m14-3.0*phi-m13);
//    		//........................................................................
//    		nn2x = ijk-strideZ*2-strideY*2;					// neighbor index (get convention)
//    		mm15 = Phi[nn2x];					// get neighbor for phi - 15
//            mm15 = 0.25*(-mm15+5.0*m15-3.0*phi-m16);
//    		//........................................................................
//    		nn2x = ijk+strideZ*2+strideY*2;					// neighbor index (get convention)
//    		mm16 = Phi[nn2x];					// get neighbor for phi - 16
//            mm16 = 0.25*(-mm16+5.0*m16-3.0*phi-m15);
//    		//........................................................................
//    		nn2x = ijk+strideZ*2-strideY*2;					// neighbor index (get convention)
//    		mm17 = Phi[nn2x];					// get neighbor for phi - 17
//            mm17 = 0.25*(-mm17+5.0*m17-3.0*phi-m18);
//    		//........................................................................
//    		nn2x = ijk-strideZ*2+strideY*2;					// neighbor index (get convention)
//    		mm18 = Phi[nn2x];					// get neighbor for phi - 18
//            mm18 = 0.25*(-mm18+5.0*m18-3.0*phi-m17);
//
//
//    		//............Compute the Color Gradient...................................
//    		nx = -3.0*1.0/18.0*(m1-m2+0.5*(m7-m8+m9-m10+m11-m12+m13-m14));
//    		ny = -3.0*1.0/18.0*(m3-m4+0.5*(m7-m8-m9+m10+m15-m16+m17-m18));
//    		nz = -3.0*1.0/18.0*(m5-m6+0.5*(m11-m12-m13+m14+m15-m16-m17+m18));
//    		//............Compute the Chemical Potential...............................
//            //chem = 2.0*3.0/18.0*(m1+m2+m3+m4+m5+m6-6*phi+0.5*(m7+m8+m9+m10+m11+m12+m13+m14+m15+m16+m17+m18-12*phi));//intermediate var, i.e. the laplacian
//            //chem = 4.0*beta*phi*(phi+1.0)*(phi-1.0)-kappa*chem;
//            chem = 2.0*3.0/18.0*(m1+m2+m3+m4+m5+m6-6*phi_temp+0.5*(m7+m8+m9+m10+m11+m12+m13+m14+m15+m16+m17+m18-12*phi_temp));//intermediate var, i.e. the laplacian
//            chem = 4.0*beta*phi_temp*(phi_temp+1.0)*(phi_temp-1.0)-kappa*chem;
//    		//............Compute the Mixed Gradient...................................
//    		mgx = -3.0*1.0/18.0*(mm1-mm2+0.5*(mm7-mm8+mm9-mm10+mm11-mm12+mm13-mm14));
//    		mgy = -3.0*1.0/18.0*(mm3-mm4+0.5*(mm7-mm8-mm9+mm10+mm15-mm16+mm17-mm18));
//    		mgz = -3.0*1.0/18.0*(mm5-mm6+0.5*(mm11-mm12-mm13+mm14+mm15-mm16-mm17+mm18));
//
//            //de-noise color gradient and mixed gradient
//            C = sqrt(nx*nx+ny*ny+nz*nz);
//            if (C<1.0e-12) nx=ny=nz=0.0;
//            double mg_mag = sqrt(mgx*mgx+mgy*mgy+mgz*mgz);
//            if (mg_mag<1.0e-12) mgx=mgy=mgz=0.0;
//            //maybe you can also de-noise chemical potential ? within the bulk phase chem should be ZERO
//            if (fabs(chem)<1.0e-12) chem=0.0;
//    		
//    		// q=0
//    		m0 = dist[n];
//    		// q=1
//    		m1 = dist[2*Np+n]; 
//
//            // q=2
//    		m2 = dist[1*Np+n];  
//
//    		// q=3
//    		m3 = dist[4*Np+n];
//
//    		// q = 4
//    		m4 = dist[3*Np+n];
//
//    		// q=5
//    		m5 = dist[6*Np+n];
//
//    		// q = 6
//    		m6 = dist[5*Np+n];
//    		
//    		// q=7
//    		m7 = dist[8*Np+n];
//
//    		// q = 8
//    		m8 = dist[7*Np+n];
//
//    		// q=9
//    		m9 = dist[10*Np+n];
//
//    		// q = 10
//    		m10 = dist[9*Np+n];
//
//    		// q=11
//    		m11 = dist[12*Np+n];
//
//    		// q=12
//    		m12 = dist[11*Np+n];
//
//    		// q=13
//    		m13 = dist[14*Np+n];
//
//    		// q=14
//    		m14 = dist[13*Np+n];
//
//    		// q=15
//    		m15 = dist[16*Np+n];
//
//    		// q=16
//    		m16 = dist[15*Np+n];
//
//    		// q=17
//    		m17 = dist[18*Np+n];
//
//    		// q=18
//    		m18 = dist[17*Np+n];
//
//            //compute fluid velocity
//            ux = 3.0/rho0*(m1-m2+m7-m8+m9-m10+m11-m12+m13-m14+0.5*(chem*nx+Fx)/3.0);
//            uy = 3.0/rho0*(m3-m4+m7-m8-m9+m10+m15-m16+m17-m18+0.5*(chem*ny+Fy)/3.0);
//            uz = 3.0/rho0*(m5-m6+m11-m12-m13+m14+m15-m16-m17+m18+0.5*(chem*nz+Fz)/3.0);
//
//            //compute pressure
//            p = (m0+m2+m1+m4+m3+m6+m5+m8+m7+m10+m9+m12+m11+m14+m13+m16+m15+m18+m17)
//                      +0.5*(rhoA-rhoB)/2.0/3.0*(ux*nx+uy*ny+uz*nz);
//
//            //compute equilibrium distributions
//                feq0 = 0.3333333333333333*p - 0.25*(Fx*ux + Fy*uy + Fz*uz)*(-0.6666666666666666 + ux*ux + uy*uy + uz*uz) - 
//         0.16666666666666666*rho0*(ux*ux + uy*uy + uz*uz) - 0.5*(-(nx*ux) - ny*uy - nz*uz)*
//          (-0.08333333333333333*(rhoA - rhoB)*(ux*ux + uy*uy + uz*uz) + chem*(0.3333333333333333 - 0.5*(ux*ux + uy*uy + uz*uz)));
//                feq1 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-ux*ux + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*uz)) - 
//         0.125*(Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 
//           0.3333333333333333*(-2.*ux + ux*ux + uy*uy + uz*uz)) - 0.0625*(nx - nx*ux - ny*uy - nz*uz)*
//          (2*chem*ux*ux - 0.3333333333333333*((-rhoA + rhoB)*ux*ux + 2*chem*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
//           0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*uz)));
//                feq2 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-ux*ux + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*uz)) - 
//         0.125*(Fx + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 
//           0.3333333333333333*(2.*ux + ux*ux + uy*uy + uz*uz)) - 0.0625*(nx + nx*ux + ny*uy + nz*uz)*
//          (-2.*chem*ux*ux + 0.1111111111111111*(-4.*chem + rhoB*(-2.*ux - 1.*ux*ux - 1.*uy*uy - 1.*uz*uz) + 
//             rhoA*(2.*ux + ux*ux + uy*uy + uz*uz)) + 0.3333333333333333*((-1.*rhoA + rhoB)*ux*ux + 
//             chem*(4.*ux + 2.*ux*ux + 2.*uy*uy + 2.*uz*uz)));
//                feq3 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uy*uy + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*uz)) - 
//         0.125*(Fx*ux + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy + 
//           0.3333333333333333*(ux*ux - 2.*uy + uy*uy + uz*uz)) - 0.0625*(ny - nx*ux - ny*uy - nz*uz)*
//          (2*chem*uy*uy - 0.3333333333333333*((-rhoA + rhoB)*uy*uy + 2*chem*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
//           0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*uz)));
//                feq4 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uy*uy + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*uz)) - 
//         0.125*(Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy + 
//           0.3333333333333333*(ux*ux + 2.*uy + uy*uy + uz*uz)) - 0.0625*(ny + nx*ux + ny*uy + nz*uz)*
//          (-2.*chem*uy*uy + 0.1111111111111111*(-4.*chem + rhoB*(-1.*ux*ux - 2.*uy - 1.*uy*uy - 1.*uz*uz) + 
//             rhoA*(ux*ux + 2.*uy + uy*uy + uz*uz)) + 0.3333333333333333*((-1.*rhoA + rhoB)*uy*uy + 
//             chem*(2.*ux*ux + 4.*uy + 2.*uy*uy + 2.*uz*uz))); 
//                feq5 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uz*uz + 0.3333333333333333*(ux*ux + uy*uy + (-2 + uz)*uz)) - 
//         0.125*(Fx*ux + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*uz*uz + 
//           0.3333333333333333*(ux*ux + uy*uy + (-2. + uz)*uz)) - 0.0625*(nx*ux + ny*uy + nz*(-1. + uz))*
//          (-2.*chem*uz*uz + 0.1111111111111111*(-4.*chem + rhoB*(-1.*ux*ux - 1.*uy*uy + (2. - 1.*uz)*uz) + 
//             rhoA*(ux*ux + uy*uy + (-2. + uz)*uz)) + 0.3333333333333333*((-1.*rhoA + rhoB)*uz*uz + 
//             chem*(2.*ux*ux + 2.*uy*uy + uz*(-4. + 2.*uz)))); 
//                feq6 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uz*uz + 0.3333333333333333*(ux*ux + uy*uy + uz*(2 + uz))) - 
//         0.125*(Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*uz*uz + 
//           0.3333333333333333*(ux*ux + uy*uy + uz*(2. + uz))) - 0.0625*(nz + nx*ux + ny*uy + nz*uz)*
//          (-2.*chem*uz*uz + 0.1111111111111111*(-4.*chem + rhoB*(-1.*ux*ux - 1.*uy*uy + (-2. - 1.*uz)*uz) + 
//             rhoA*(ux*ux + uy*uy + uz*(2. + uz))) + 0.3333333333333333*((-1.*rhoA + rhoB)*uz*uz + 
//             chem*(2.*ux*ux + 2.*uy*uy + uz*(4. + 2.*uz)))); 
//                feq7 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(ux + uy)*(ux + uy) + 0.3333333333333333*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) - 
//         0.0625*(Fx*(-1. + ux) + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uy - 1.*uy*uy + 
//           0.3333333333333333*(-2.*ux + ux*ux - 2.*uy + uy*uy + uz*uz)) - 0.03125*(nx + ny - nx*ux - ny*uy - nz*uz)*
//          (2*chem*(ux + uy)*(ux + uy) + 0.3333333333333333*((rhoA - rhoB)*(ux + uy)*(ux + uy) - 2*chem*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
//           0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz))); 
//                feq8 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(ux + uy)*(ux + uy) + 0.3333333333333333*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) - 
//         0.0625*(Fx + Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uy - 1.*uy*uy + 
//           0.3333333333333333*(2.*ux + ux*ux + 2.*uy + uy*uy + uz*uz)) - 0.03125*(-(nx*(1 + ux)) - ny*(1 + uy) - nz*uz)*
//          (2*chem*(ux + uy)*(ux + uy) - 0.3333333333333333*(-((rhoA - rhoB)*(ux + uy)*(ux + uy)) + 
//             2*chem*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 0.1111111111111111*
//            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz))); 
//                feq9 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(ux - uy)*(ux - uy) + 0.3333333333333333*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) - 
//         0.0625*(Fy + Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uy - 1.*uy*uy + 
//           0.3333333333333333*(-2.*ux + ux*ux + 2.*uy + uy*uy + uz*uz)) - 0.03125*(nx - nx*ux - ny*(1 + uy) - nz*uz)*
//          (2*chem*(ux - uy)*(ux - uy) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uy)*(ux - uy)) + 
//             2*chem*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 0.1111111111111111*
//            (4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz))); 
//                feq10 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(ux - uy)*(ux - uy) + 0.3333333333333333*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) - 
//         0.0625*(Fx*(1 + ux) + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uy - 1.*uy*uy + 
//           0.3333333333333333*(2.*ux + ux*ux - 2.*uy + uy*uy + uz*uz)) - 0.03125*(ny - nx*(1 + ux) - ny*uy - nz*uz)*
//          (2*chem*(ux - uy)*(ux - uy) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uy)*(ux - uy)) + 
//             2*chem*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 0.1111111111111111*
//            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz))); 
//                feq11 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(ux + uz)*(ux + uz) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) - 
//         0.0625*(Fx*(-1. + ux) + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uz - 1.*uz*uz + 
//           0.3333333333333333*(-2.*ux + ux*ux + uy*uy + (-2. + uz)*uz)) - 0.03125*(nx + nz - nx*ux - ny*uy - nz*uz)*
//          (2*chem*(ux + uz)*(ux + uz) + 0.3333333333333333*((rhoA - rhoB)*(ux + uz)*(ux + uz) - 2*chem*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
//           0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz))); 
//                feq12 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(ux + uz)*(ux + uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) - 
//         0.0625*(Fx + Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uz - 1.*uz*uz + 
//           0.3333333333333333*(2.*ux + ux*ux + uy*uy + uz*(2. + uz))) - 0.03125*(-(nx*(1 + ux)) - ny*uy - nz*(1 + uz))*
//          (2*chem*(ux + uz)*(ux + uz) - 0.3333333333333333*(-((rhoA - rhoB)*(ux + uz)*(ux + uz)) + 
//             2*chem*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
//            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))); 
//                feq13 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(ux - uz)*(ux - uz) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) - 
//         0.0625*(Fz + Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uz - 1.*uz*uz + 
//           0.3333333333333333*(-2.*ux + ux*ux + uy*uy + uz*(2. + uz))) - 0.03125*(nx - nx*ux - ny*uy - nz*(1 + uz))*
//          (2*chem*(ux - uz)*(ux - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uz)*(ux - uz)) + 
//             2*chem*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
//            (4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*(2 + uz)))); 
//                feq14 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(ux - uz)*(ux - uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) - 
//         0.0625*(Fx*(1 + ux) + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uz - 1.*uz*uz + 
//           0.3333333333333333*(2.*ux + ux*ux + uy*uy + (-2. + uz)*uz)) - 0.03125*(nz - nx*(1 + ux) - ny*uy - nz*uz)*
//          (2*chem*(ux - uz)*(ux - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uz)*(ux - uz)) + 
//             2*chem*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 0.1111111111111111*
//            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz))); 
//                feq15 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) - 
//         0.0625*(Fx*ux + Fy*(-1. + uy) + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*uy*uy - 2.*uy*uz - 1.*uz*uz + 
//           0.3333333333333333*(ux*ux - 2.*uy + uy*uy + (-2. + uz)*uz)) - 0.03125*(ny + nz - nx*ux - ny*uy - nz*uz)*
//          (2*chem*(uy + uz)*(uy + uz) + 0.3333333333333333*((rhoA - rhoB)*(uy + uz)*(uy + uz) - 2*chem*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) + 
//           0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))); 
//                feq16 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) - 
//         0.0625*(Fy + Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy - 2.*uy*uz - 1.*uz*uz + 
//           0.3333333333333333*(ux*ux + 2.*uy + uy*uy + uz*(2. + uz))) - 0.03125*(-(nx*ux) - ny*(1 + uy) - nz*(1 + uz))*
//          (2*chem*(uy + uz)*(uy + uz) - 0.3333333333333333*(-((rhoA - rhoB)*(uy + uz)*(uy + uz)) + 
//             2*chem*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
//            (4*chem - (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))); 
//                feq17 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) - 
//         0.0625*(Fz + Fx*ux + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy + 2.*uy*uz - 1.*uz*uz + 
//           0.3333333333333333*(ux*ux - 2.*uy + uy*uy + uz*(2. + uz))) - 0.03125*(ny - nx*ux - ny*uy - nz*(1 + uz))*
//          (2*chem*(uy - uz)*(uy - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(uy - uz)*(uy - uz)) + 
//             2*chem*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
//            (4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))); 
//                feq18 = 0.027777777777777776*p - 0.041666666666666664*rho0*
//          (-(uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) - 
//         0.0625*(Fx*ux + Fy*(1 + uy) + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*uy*uy + 2.*uy*uz - 1.*uz*uz + 
//           0.3333333333333333*(ux*ux + 2.*uy + uy*uy + (-2. + uz)*uz)) - 0.03125*(nz - nx*ux - ny*(1 + uy) - nz*uz)*
//          (2*chem*(uy - uz)*(uy - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(uy - uz)*(uy - uz)) + 
//             2*chem*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 0.1111111111111111*
//            (4*chem - (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz))); 
//
//            //------------------------------------------------- BCK collison ------------------------------------------------------------//
//    		// q=0
//    		dist[n] = m0 - (m0-feq0)/tau + 0.25*(2*(Fx*ux + Fy*uy + Fz*uz)*(-0.6666666666666666 + ux*ux + uy*uy + uz*uz) + 
//         (mgx*ux + mgy*uy + mgz*uz)*(2*chem*(ux*ux + uy*uy + uz*uz) + 
//            0.3333333333333333*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + uz*uz))));
//
//    		// q = 1
//    		dist[1*Np+n] = m1 - (m1-feq1)/tau + 0.125*(2*(Fx*(-1 + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - ux*ux + 
//            0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
//         (mgx*(-1 + ux) + mgy*uy + mgz*uz)*(-2*chem*(ux*ux) + 
//            0.3333333333333333*((-rhoA + rhoB)*(ux*ux) + 2*chem*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*uz))));
//
//    		// q=2
//    		dist[2*Np+n] = m2 - (m2-feq2)/tau + 0.125*(2*(Fx + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - ux*ux + 
//            0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*uz)) + 
//         (mgx + mgx*ux + mgy*uy + mgz*uz)*(-2*chem*(ux*ux) + 
//            0.3333333333333333*((-rhoA + rhoB)*(ux*ux) + 2*chem*(2*ux + ux*ux + uy*uy + uz*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*uz))));
//
//    		// q = 3
//    		dist[3*Np+n] = m3 - (m3-feq3)/tau + 0.125*(2*(Fx*ux + Fy*(-1 + uy) + Fz*uz)*(-0.2222222222222222 - uy*uy + 
//            0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
//         (mgx*ux + mgy*(-1 + uy) + mgz*uz)*(-2*chem*(uy*uy) + 
//            0.3333333333333333*((-rhoA + rhoB)*(uy*uy) + 2*chem*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*uz))));
//
//    		// q = 4
//    		dist[4*Np+n] = m4 - (m4-feq4)/tau + 0.125*(2*(Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - uy*uy + 
//            0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*uz)) + 
//         (mgy + mgx*ux + mgy*uy + mgz*uz)*(-2*chem*(uy*uy) + 
//            0.3333333333333333*((-rhoA + rhoB)*(uy*uy) + 2*chem*(ux*ux + 2*uy + uy*uy + uz*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*uz))));
//
//    		// q = 5
//    		dist[5*Np+n] = m5 - (m5-feq5)/tau + 0.125*(2*(Fx*ux + Fy*uy + Fz*(-1 + uz))*(-0.2222222222222222 - uz*uz + 
//            0.3333333333333333*(ux*ux + uy*uy + (-2 + uz)*uz)) + 
//         (mgx*ux + mgy*uy + mgz*(-1 + uz))*(-2*chem*(uz*uz) + 
//            0.3333333333333333*((-rhoA + rhoB)*(uz*uz) + 2*chem*(ux*ux + uy*uy + (-2 + uz)*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + (-2 + uz)*uz))));
//
//    		// q = 6
//    		dist[6*Np+n] = m6 - (m6-feq6)/tau + 0.125*(2*(Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - uz*uz + 
//            0.3333333333333333*(ux*ux + uy*uy + uz*(2 + uz))) + 
//         (mgz + mgx*ux + mgy*uy + mgz*uz)*(-2*chem*(uz*uz) + 
//            0.3333333333333333*((-rhoA + rhoB)*(uz*uz) + 2*chem*(ux*ux + uy*uy + uz*(2 + uz))) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + uz*(2 + uz)))));
//
//    		// q = 7
//    		dist[7*Np+n] = m7 - (m7-feq7)/tau + 0.0625*(-2*(Fx*(-1 + ux) + Fy*(-1 + uy) + Fz*uz)*
//          (0.2222222222222222 + (ux + uy)*(ux + uy) - 
//            0.3333333333333333*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
//         (mgx*(-1 + ux) + mgy*(-1 + uy) + mgz*uz)*
//          (-2*chem*((ux + uy)*(ux + uy)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((ux + uy)*(ux + uy))) + 2*chem*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz))));
//
//    		// q = 8
//    		dist[8*Np+n] = m8 - (m8-feq8)/tau + 0.0625*(2*(Fx + Fy + Fx*ux + Fy*uy + Fz*uz)*
//          (-0.2222222222222222 - (ux + uy)*(ux + uy) + 
//            0.3333333333333333*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
//         (mgx + mgy + mgx*ux + mgy*uy + mgz*uz)*
//          (-2*chem*((ux + uy)*(ux + uy)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((ux + uy)*(ux + uy))) + 2*chem*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz))));
//
//    		// q = 9
//    		dist[9*Np+n] = m9 - (m9-feq9)/tau + 0.0625*(2*(Fy + Fx*(-1 + ux) + Fy*uy + Fz*uz)*
//          (-0.2222222222222222 - (ux - uy)*(ux - uy) + 
//            0.3333333333333333*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
//         (mgy + mgx*(-1 + ux) + mgy*uy + mgz*uz)*
//          (-2*chem*((ux - uy)*(ux - uy)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((ux - uy)*(ux - uy))) + 2*chem*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz))));
//
//    		// q = 10
//    		dist[10*Np+n] = m10 - (m10-feq10)/tau + 0.0625*(2*(Fx*(1 + ux) + Fy*(-1 + uy) + Fz*uz)*
//          (-0.2222222222222222 - (ux - uy)*(ux - uy) + 
//            0.3333333333333333*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
//         (mgx*(1 + ux) + mgy*(-1 + uy) + mgz*uz)*
//          (-2*chem*((ux - uy)*(ux - uy)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((ux - uy)*(ux - uy))) + 2*chem*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz))));
//
//    		// q = 11
//    		dist[11*Np+n] = m11 - (m11-feq11)/tau + 0.0625*(-2*(Fx*(-1 + ux) + Fy*uy + Fz*(-1 + uz))*
//          (0.2222222222222222 + (ux + uz)*(ux + uz) - 
//            0.3333333333333333*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
//         (mgx*(-1 + ux) + mgy*uy + mgz*(-1 + uz))*
//          (-2*chem*((ux + uz)*(ux + uz)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((ux + uz)*(ux + uz))) + 2*chem*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz))));
//
//    		// q = 12
//    		dist[12*Np+n] = m12 - (m12-feq12)/tau + 0.0625*(2*(Fx + Fz + Fx*ux + Fy*uy + Fz*uz)*
//          (-0.2222222222222222 - (ux + uz)*(ux + uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))
//           + (mgx + mgz + mgx*ux + mgy*uy + mgz*uz)*
//          (-2*chem*((ux + uz)*(ux + uz)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((ux + uz)*(ux + uz))) + 2*chem*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))));
//
//    		// q = 13
//    		dist[13*Np+n] = m13 - (m13-feq13)/tau + 0.0625*(2*(Fz + Fx*(-1 + ux) + Fy*uy + Fz*uz)*
//          (-0.2222222222222222 - (ux - uz)*(ux - uz) + 
//            0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
//         (mgz + mgx*(-1 + ux) + mgy*uy + mgz*uz)*
//          (-2*chem*((ux - uz)*(ux - uz)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((ux - uz)*(ux - uz))) + 2*chem*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*(2 + uz)))));
//
//    		// q= 14
//    		dist[14*Np+n] = m14 - (m14-feq14)/tau + 0.0625*(2*(Fx*(1 + ux) + Fy*uy + Fz*(-1 + uz))*
//          (-0.2222222222222222 - (ux - uz)*(ux - uz) + 
//            0.3333333333333333*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
//         (mgx*(1 + ux) + mgy*uy + mgz*(-1 + uz))*
//          (-2*chem*((ux - uz)*(ux - uz)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((ux - uz)*(ux - uz))) + 2*chem*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz))));
//
//    		// q = 15
//    		dist[15*Np+n] = m15 - (m15-feq15)/tau + 0.0625*(-2*(Fx*ux + Fy*(-1 + uy) + Fz*(-1 + uz))*
//          (0.2222222222222222 + (uy + uz)*(uy + uz) - 0.3333333333333333*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))
//           + (mgx*ux + mgy*(-1 + uy) + mgz*(-1 + uz))*
//          (-2*chem*((uy + uz)*(uy + uz)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((uy + uz)*(uy + uz))) + 2*chem*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))));
//
//    		// q = 16
//    		dist[16*Np+n] = m16 - (m16-feq16)/tau + 0.0625*(2*(Fy + Fz + Fx*ux + Fy*uy + Fz*uz)*
//          (-0.2222222222222222 - (uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))
//           + (mgy + mgz + mgx*ux + mgy*uy + mgz*uz)*
//          (-2*chem*((uy + uz)*(uy + uz)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((uy + uz)*(uy + uz))) + 2*chem*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))));
//
//    		// q = 17
//    		dist[17*Np+n] = m17 - (m17-feq17)/tau + 0.0625*(2*(Fz + Fx*ux + Fy*(-1 + uy) + Fz*uz)*
//          (-0.2222222222222222 - (uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))
//           + (mgz + mgx*ux + mgy*(-1 + uy) + mgz*uz)*
//          (-2*chem*((uy - uz)*(uy - uz)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((uy - uz)*(uy - uz))) + 2*chem*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))));
//
//    		// q = 18
//    		dist[18*Np+n] = m18 - (m18-feq18)/tau + 0.0625*(2*(Fx*ux + Fy*(1 + uy) + Fz*(-1 + uz))*
//          (-0.2222222222222222 - (uy - uz)*(uy - uz) + 
//            0.3333333333333333*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 
//         (mgx*ux + mgy*(1 + uy) + mgz*(-1 + uz))*
//          (-2*chem*((uy - uz)*(uy - uz)) + 0.3333333333333333*
//             (-((rhoA - rhoB)*((uy - uz)*(uy - uz))) + 2*chem*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 
//            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz))));
//            //----------------------------------------------------------------------------------------------------------------------------------------//
//
//            // ----------------------------- compute phase field evolution ----------------------------------------
//            //Normalize the Color Gradient
//            C = sqrt(nx*nx+ny*ny+nz*nz);
//            double ColorMag = C;
//            if (C==0.0) ColorMag=1.0;
//            nx = nx/ColorMag;
//            ny = ny/ColorMag;
//            nz = nz/ColorMag;		
//            //compute surface tension-related parameter
//            //theta = 4.5*M*2.0*(1-phi*phi)/W;
//            theta = 4.5*M*2.0*(1-phi_temp*phi_temp)/W;
//
//            //load distributions of phase field
//            //q=0
//            h0 = hq[n];
//            //q=1
//            h1 = hq[2*Np+n]; 
//
//            //q=2
//            h2 = hq[1*Np+n];  
//
//            //q=3
//            h3 = hq[4*Np+n];
//
//            //q=4
//            h4 = hq[3*Np+n];
//
//            //q=5
//            h5 = hq[6*Np+n];
//
//            //q=6
//            h6 = hq[5*Np+n];
//
//            //-------------------------------- BGK collison for phase field ---------------------------------//
//            // q = 0
//            hq[n] = h0 - (h0 - 0.3333333333333333*phi)/tauM;
//
//            // q = 1
//            hq[1*Np+n] = h1 - (h1 - 0.1111111111111111*nx*theta - phi*(0.1111111111111111 + 0.5*ux))/tauM;
//
//            // q = 2
//            hq[2*Np+n] = h2 - (h2 + 0.1111111111111111*nx*theta - phi*(0.1111111111111111 - 0.5*ux))/tauM;
//
//            // q = 3
//            hq[3*Np+n] = h3 - (h3 - 0.1111111111111111*ny*theta - phi*(0.1111111111111111 + 0.5*uy))/tauM;
//
//            // q = 4
//            hq[4*Np+n] = h4 - (h4 + 0.1111111111111111*ny*theta - phi*(0.1111111111111111 - 0.5*uy))/tauM;
//
//            // q = 5
//            hq[5*Np+n] = h5 - (h5 - 0.1111111111111111*nz*theta - phi*(0.1111111111111111 + 0.5*uz))/tauM;
//
//            // q = 6
//            hq[6*Np+n] = h6 - (h6 + 0.1111111111111111*nz*theta - phi*(0.1111111111111111 - 0.5*uz))/tauM;
//            //........................................................................
//
//            //Update velocity on device
//    		Vel[0*Np+n] = ux;
//    		Vel[1*Np+n] = uy;
//    		Vel[2*Np+n] = uz;
//            //Update pressure on device
//            Pressure[n] = p;
//            //Update chemical potential on device
//            mu_phi[n] = chem;
//            //Update color gradient on device
//    		ColorGrad[0*Np+n] = nx;
//    		ColorGrad[1*Np+n] = ny;
//    		ColorGrad[2*Np+n] = nz;
//
//    	}
//	}
//}

__global__ void dvc_ScaLBL_D3Q19_AAodd_FreeLeeModel_Combined_HigherOrder(int *neighborList, int *Map, double *dist, double *hq, double *Den,	double *Phi, double *mu_phi, double *Vel, double *Pressure, double *ColorGrad, 
        double rhoA, double rhoB, double tauA, double tauB, double tauM, double kappa, double beta, double W, double Fx, double Fy, double Fz, 
        int strideY, int strideZ, int start, int finish, int Np){

	int n,nn,ijk;
	int nr1,nr2,nr3,nr4,nr5,nr6,nr7,nr8,nr9,nr10,nr11,nr12,nr13,nr14,nr15,nr16,nr17,nr18;
    double ux,uy,uz;//fluid velocity 
    double p;//pressure
    double chem;//chemical potential
    double phi; //phase field
    double rho0;//fluid density
	// register variables for neighbors of phase field
    double m0;
    // w(|c|^2) = 1 
	double m1,m2,m3,m4,m5,m6;
    // w(|c|^2) = 2
	double m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
    // w(|c|^2) = 3
    double m19,m20,m21,m22,m23,m24,m25,m26;
    // w(|c|^2) = 4
    double m27,m28,m29,m30,m31,m32;
    // w(|c|^2) = 5
    double m33,m34,m35,m36,m37,m38,m39,m40,m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,m51,m52,m53,m54,m55,m56;
    // w(|c|^2) = 6
    double m57,m58,m59,m60,m61,m62,m63,m64,m65,m66,m67,m68,m69,m70,m71,m72,m73,m74,m75,m76,m77,m78,m79,m80;
    // w(|c|^2) = 8
    double m81,m82,m83,m84,m85,m86,m87,m88,m89,m90,m91,m92;
	//double mm1,mm2,mm4,mm6,mm8,mm9,mm10,mm11,mm12,mm13,mm14,mm15,mm16,mm17,mm18;
	//double mm3,mm5,mm7;
    double feq0,feq1,feq2,feq3,feq4,feq5,feq6,feq7,feq8,feq9,feq10,feq11,feq12,feq13,feq14,feq15,feq16,feq17,feq18;
    double nx,ny,nz;//normal color gradient
    //double mgx,mgy,mgz;//mixed gradient reaching secondary neighbor

    double h0,h1,h2,h3,h4,h5,h6;//distributions for LB phase field
	double tau;//position dependent LB relaxation time for fluid
    double C,theta;
     double M = 2.0/9.0*(tauM-0.5);//diffusivity (or mobility) for the phase field D3Q7
    double phi_temp;

    int S = Np/NBLOCKS/NTHREADS + 1;
    for (int s=0; s<S; s++){
    	//........Get 1-D index for this thread....................
    	n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

    	if ( n<finish ){
    		rho0 = Den[n];//load density

    		// Get the 1D index based on regular data layout
    		ijk = Map[n];
            phi = Phi[ijk];// load phase field
            phi_temp = phi;
            if (phi>1.f) phi_temp=1.0;
            if (phi<-1.f) phi_temp=-1.0;

    		// local relaxation time
    		tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);

    		//					COMPUTE THE COLOR GRADIENT
    		//........................................................................
    		//.................Read Phase Indicator Values............................
    		//........................................................................
            //-------------------
            //--- w(|c|^2) = 1--- 
            //-------------------
    		nn = ijk-1;							// neighbor index (get convention)
    		m2 = Phi[nn];						// get neighbor for phi - 1
    		//........................................................................
    		nn = ijk+1;							// neighbor index (get convention)
    		m1 = Phi[nn];						// get neighbor for phi - 2
    		//........................................................................
    		nn = ijk-strideY;							// neighbor index (get convention)
    		m4 = Phi[nn];					// get neighbor for phi - 3
    		//........................................................................
    		nn = ijk+strideY;							// neighbor index (get convention)
    		m3 = Phi[nn];					// get neighbor for phi - 4
    		//........................................................................
    		nn = ijk-strideZ;						// neighbor index (get convention)
    		m6 = Phi[nn];					// get neighbor for phi - 5
    		//........................................................................
    		nn = ijk+strideZ;						// neighbor index (get convention)
    		m5 = Phi[nn];					// get neighbor for phi - 6
    		//........................................................................

            //-------------------
            //--- w(|c|^2) = 2--- 
            //-------------------
    		nn = ijk-strideY-1;						// neighbor index (get convention)
    		m8 = Phi[nn];					// get neighbor for phi - 7
    		//........................................................................
    		nn = ijk+strideY+1;						// neighbor index (get convention)
    		m7 = Phi[nn];					// get neighbor for phi - 8
    		//........................................................................
    		nn = ijk+strideY-1;						// neighbor index (get convention)
    		m10 = Phi[nn];					// get neighbor for phi - 9
    		//........................................................................
    		nn = ijk-strideY+1;						// neighbor index (get convention)
    		m9 = Phi[nn];					// get neighbor for phi - 10
    		//........................................................................
    		nn = ijk-strideZ-1;						// neighbor index (get convention)
    		m12 = Phi[nn];					// get neighbor for phi - 11
    		//........................................................................
    		nn = ijk+strideZ+1;						// neighbor index (get convention)
    		m11 = Phi[nn];					// get neighbor for phi - 12
    		//........................................................................
    		nn = ijk+strideZ-1;						// neighbor index (get convention)
    		m14 = Phi[nn];					// get neighbor for phi - 13
    		//........................................................................
    		nn = ijk-strideZ+1;						// neighbor index (get convention)
    		m13 = Phi[nn];					// get neighbor for phi - 14
    		//........................................................................
    		nn = ijk-strideZ-strideY;					// neighbor index (get convention)
    		m16 = Phi[nn];					// get neighbor for phi - 15
    		//........................................................................
    		nn = ijk+strideZ+strideY;					// neighbor index (get convention)
    		m15 = Phi[nn];					// get neighbor for phi - 16
    		//........................................................................
    		nn = ijk+strideZ-strideY;					// neighbor index (get convention)
    		m18 = Phi[nn];					// get neighbor for phi - 17
    		//........................................................................
    		nn = ijk-strideZ+strideY;					// neighbor index (get convention)
    		m17 = Phi[nn];					// get neighbor for phi - 18
    		//........................................................................

            //-------------------
            //--- w(|c|^2) = 3--- 
            //-------------------
            // c19 = (1,1,1) = (cx,cy,cz)
    		nn = ijk+strideZ+strideY+1;	
    		m19 = Phi[nn];			
    		//........................................................................
            // c20 = (-1,-1,-1)
    		nn = ijk-strideZ-strideY-1;	
    		m20 = Phi[nn];			
    		//........................................................................
            // c21 = (1,-1,1)
    		nn = ijk+strideZ-strideY+1;	
    		m21 = Phi[nn];			
    		//........................................................................
            // c22 = (-1,1,-1)
    		nn = ijk-strideZ+strideY-1;	
    		m22 = Phi[nn];			
    		//........................................................................
            // c23 = (1,1,-1)
    		nn = ijk-strideZ+strideY+1;	
    		m23 = Phi[nn];			
    		//........................................................................
            // c24 = (-1,-1,1)
    		nn = ijk+strideZ-strideY-1;	
    		m24 = Phi[nn];			
    		//........................................................................
            // c25 = (-1,1,1)
    		nn = ijk+strideZ+strideY-1;	
    		m25 = Phi[nn];			
    		//........................................................................
            // c26 = (1,-1,-1)
    		nn = ijk-strideZ-strideY+1;	
    		m26 = Phi[nn];			
    		//........................................................................
            
            //-------------------
            //--- w(|c|^2) = 4--- 
            //-------------------
            // c27 = (2,0,0) = (cx,cy,cz)
    		nn = ijk+2;		
    		m27 = Phi[nn];			
    		//........................................................................
            // c28 = (-2,0,0) = (cx,cy,cz)
    		nn = ijk-2;		
    		m28 = Phi[nn];			
    		//........................................................................
            // c29 = (0,2,0) = (cx,cy,cz)
    		nn = ijk+2*strideY;		
    		m29 = Phi[nn];			
    		//........................................................................
            // c30 = (0,-2,0) = (cx,cy,cz)
    		nn = ijk-2*strideY;		
    		m30 = Phi[nn];			
    		//........................................................................
            // c31 = (0,0,2) = (cx,cy,cz)
    		nn = ijk+2*strideZ;		
    		m31 = Phi[nn];			
    		//........................................................................
            // c32 = (0,0,-2) = (cx,cy,cz)
    		nn = ijk-2*strideZ;		
    		m32 = Phi[nn];			
    		//........................................................................

            //-------------------
            //--- w(|c|^2) = 5--- 
            //-------------------
            // c33 = (1,2,0) = (cx,cy,cz)
    		nn = ijk+2*strideY+1;				
    		m33 = Phi[nn];					
    		//........................................................................
            // c34 = (-1,-2,0) = (cx,cy,cz)
    		nn = ijk-2*strideY-1;				
    		m34 = Phi[nn];					
    		//........................................................................
            // c35 = (2,1,0) = (cx,cy,cz)
    		nn = ijk+strideY+2;				
    		m35 = Phi[nn];					
    		//........................................................................
            // c36 = (-2,-1,0) = (cx,cy,cz)
    		nn = ijk-strideY-2;				
    		m36 = Phi[nn];					
    		//........................................................................
            // c37 = (1,-2,0) = (cx,cy,cz)
    		nn = ijk-2*strideY+1;				
    		m37 = Phi[nn];					
    		//........................................................................
            // c38 = (-1,2,0) = (cx,cy,cz)
    		nn = ijk+2*strideY-1;				
    		m38 = Phi[nn];					
    		//........................................................................
            // c39 = (2,-1,0) = (cx,cy,cz)
    		nn = ijk-strideY+2;				
    		m39 = Phi[nn];					
    		//........................................................................
            // c40 = (-2,1,0) = (cx,cy,cz)
    		nn = ijk+strideY-2;				
    		m40 = Phi[nn];					
    		//........................................................................
            // c41 = (1,0,2) = (cx,cy,cz)
    		nn = ijk+2*strideZ+1;				
    		m41 = Phi[nn];					
    		//........................................................................
            // c42 = (-1,0,-2) = (cx,cy,cz)
    		nn = ijk-2*strideZ-1;				
    		m42 = Phi[nn];					
    		//........................................................................
            // c43 = (2,0,1) = (cx,cy,cz)
    		nn = ijk+strideZ+2;				
    		m43 = Phi[nn];					
    		//........................................................................
            // c44 = (-2,0,-1) = (cx,cy,cz)
    		nn = ijk-strideZ-2;				
    		m44 = Phi[nn];					
    		//........................................................................
            // c45 = (1,0,-2) = (cx,cy,cz)
    		nn = ijk-2*strideZ+1;				
    		m45 = Phi[nn];					
    		//........................................................................
            // c46 = (-1,0,2) = (cx,cy,cz)
    		nn = ijk+2*strideZ-1;				
    		m46 = Phi[nn];					
    		//........................................................................
            // c47 = (2,0,-1) = (cx,cy,cz)
    		nn = ijk-strideZ+2;				
    		m47 = Phi[nn];					
    		//........................................................................
            // c48 = (-2,0,1) = (cx,cy,cz)
    		nn = ijk+strideZ-2;				
    		m48 = Phi[nn];					
    		//........................................................................
            // c49 = (0,1,2) = (cx,cy,cz)
    		nn = ijk+2*strideZ+strideY;				
    		m49 = Phi[nn];					
    		//........................................................................
            // c50 = (0,-1,-2) = (cx,cy,cz)
    		nn = ijk-2*strideZ-strideY;				
    		m50 = Phi[nn];					
    		//........................................................................
            // c51 = (0,2,1) = (cx,cy,cz)
    		nn = ijk+strideZ+2*strideY;				
    		m51 = Phi[nn];					
    		//........................................................................
            // c52 = (0,-2,-1) = (cx,cy,cz)
    		nn = ijk-strideZ-2*strideY;				
    		m52 = Phi[nn];					
    		//........................................................................
            // c53 = (0,1,-2) = (cx,cy,cz)
    		nn = ijk-2*strideZ+strideY;				
    		m53 = Phi[nn];					
    		//........................................................................
            // c54 = (0,-1,2) = (cx,cy,cz)
    		nn = ijk+2*strideZ-strideY;				
    		m54 = Phi[nn];					
    		//........................................................................
            // c55 = (0,2,-1) = (cx,cy,cz)
    		nn = ijk-strideZ+2*strideY;				
    		m55 = Phi[nn];					
    		//........................................................................
            // c56 = (0,-2,1) = (cx,cy,cz)
    		nn = ijk+strideZ-2*strideY;				
    		m56 = Phi[nn];					
    		//........................................................................

            //-------------------
            //--- w(|c|^2) = 6--- 
            //-------------------
            // c57 = (1,1,2) = (cx,cy,cz)
    		nn = ijk+1+strideY+2*strideZ;				
    		m57 = Phi[nn];					
    		//........................................................................
            // c58 = (-1,-1,-2) = (cx,cy,cz)
    		nn = ijk-1-strideY-2*strideZ;				
    		m58 = Phi[nn];					
    		//........................................................................
            // c59 = (1,2,1) = (cx,cy,cz)
    		nn = ijk+1+2*strideY+strideZ;				
    		m59 = Phi[nn];					
    		//........................................................................
            // c60 = (-1,-2,-1) = (cx,cy,cz)
    		nn = ijk-1-2*strideY-strideZ;				
    		m60 = Phi[nn];					
    		//........................................................................
            // c61 = (2,1,1) = (cx,cy,cz)
    		nn = ijk+2+strideY+strideZ;				
    		m61 = Phi[nn];					
    		//........................................................................
            // c62 = (-2,-1,-1) = (cx,cy,cz)
    		nn = ijk-2-strideY-strideZ;				
    		m62 = Phi[nn];					
    		//........................................................................
            // c63 = (-1,1,2) = (cx,cy,cz)
    		nn = ijk-1+strideY+2*strideZ;				
    		m63 = Phi[nn];					
    		//........................................................................
            // c64 = (1,-1,-2) = (cx,cy,cz)
    		nn = ijk+1-strideY-2*strideZ;				
    		m64 = Phi[nn];					
    		//........................................................................
            // c65 = (1,-1,2) = (cx,cy,cz)
    		nn = ijk+1-strideY+2*strideZ;				
    		m65 = Phi[nn];					
    		//........................................................................
            // c66 = (-1,1,-2) = (cx,cy,cz)
    		nn = ijk-1+strideY-2*strideZ;				
    		m66 = Phi[nn];					
    		//........................................................................
            // c67 = (1,1,-2) = (cx,cy,cz)
    		nn = ijk+1+strideY-2*strideZ;				
    		m67 = Phi[nn];					
    		//........................................................................
            // c68 = (-1,-1,2) = (cx,cy,cz)
    		nn = ijk-1-strideY+2*strideZ;				
    		m68 = Phi[nn];					
    		//........................................................................
            // c69 = (-1,2,1) = (cx,cy,cz)
    		nn = ijk-1+2*strideY+strideZ;				
    		m69 = Phi[nn];					
    		//........................................................................
            // c70 = (1,-2,-1) = (cx,cy,cz)
    		nn = ijk+1-2*strideY-strideZ;				
    		m70 = Phi[nn];					
    		//........................................................................
            // c71 = (1,-2,1) = (cx,cy,cz)
    		nn = ijk+1-2*strideY+strideZ;				
    		m71 = Phi[nn];					
    		//........................................................................
            // c72 = (-1,2,-1) = (cx,cy,cz)
    		nn = ijk-1+2*strideY-strideZ;				
    		m72 = Phi[nn];					
    		//........................................................................
            // c73 = (1,2,-1) = (cx,cy,cz)
    		nn = ijk+1+2*strideY-strideZ;				
    		m73 = Phi[nn];					
    		//........................................................................
            // c74 = (-1,-2,1) = (cx,cy,cz)
    		nn = ijk-1-2*strideY+strideZ;				
    		m74 = Phi[nn];					
    		//........................................................................
            // c75 = (-2,1,1) = (cx,cy,cz)
    		nn = ijk-2+strideY+strideZ;				
    		m75 = Phi[nn];					
    		//........................................................................
            // c76 = (2,-1,-1) = (cx,cy,cz)
    		nn = ijk+2-strideY-strideZ;				
    		m76 = Phi[nn];					
    		//........................................................................
            // c77 = (2,-1,1) = (cx,cy,cz)
    		nn = ijk+2-strideY+strideZ;				
    		m77 = Phi[nn];					
    		//........................................................................
            // c78 = (-2,1,-1) = (cx,cy,cz)
    		nn = ijk-2+strideY-strideZ;				
    		m78 = Phi[nn];					
    		//........................................................................
            // c79 = (2,1,-1) = (cx,cy,cz)
    		nn = ijk+2+strideY-strideZ;				
    		m79 = Phi[nn];					
    		//........................................................................
            // c80 = (-2,-1,1) = (cx,cy,cz)
    		nn = ijk-2-strideY+strideZ;				
    		m80 = Phi[nn];					
    		//........................................................................

            //-------------------
            //--- w(|c|^2) = 8--- 
            //-------------------
            // c81 = (2,2,0) = (cx,cy,cz)
    		nn = ijk+2+2*strideY;				
    		m81 = Phi[nn];					
    		//........................................................................
            // c82 = (-2,-2,0) = (cx,cy,cz)
    		nn = ijk-2-2*strideY;				
    		m82 = Phi[nn];					
    		//........................................................................
            // c83 = (2,-2,0) = (cx,cy,cz)
    		nn = ijk+2-2*strideY;				
    		m83 = Phi[nn];					
    		//........................................................................
            // c84 = (-2,2,0) = (cx,cy,cz)
    		nn = ijk-2+2*strideY;				
    		m84 = Phi[nn];					
    		//........................................................................
            // c85 = (2,0,2) = (cx,cy,cz)
    		nn = ijk+2+2*strideZ;				
    		m85 = Phi[nn];					
    		//........................................................................
            // c86 = (-2,0,-2) = (cx,cy,cz)
    		nn = ijk-2-2*strideZ;				
    		m86 = Phi[nn];					
    		//........................................................................
            // c87 = (2,0,-2) = (cx,cy,cz)
    		nn = ijk+2-2*strideZ;				
    		m87 = Phi[nn];					
    		//........................................................................
            // c88 = (-2,0,2) = (cx,cy,cz)
    		nn = ijk-2+2*strideZ;				
    		m88 = Phi[nn];					
    		//........................................................................
            // c89 = (0,2,2) = (cx,cy,cz)
    		nn = ijk+2*strideY+2*strideZ;				
    		m89 = Phi[nn];					
    		//........................................................................
            // c90 = (0,-2,-2) = (cx,cy,cz)
    		nn = ijk-2*strideY-2*strideZ;				
    		m90 = Phi[nn];					
    		//........................................................................
            // c91 = (0,2,-2) = (cx,cy,cz)
    		nn = ijk+2*strideY-2*strideZ;				
    		m91 = Phi[nn];					
    		//........................................................................
            // c92 = (0,-2,2) = (cx,cy,cz)
    		nn = ijk-2*strideY+2*strideZ;				
    		m92 = Phi[nn];					
    		//........................................................................

    		//............Compute the Color Gradient...................................
            nx =  4.0/45.0*(m1-m2)
                + 1.0/21.0*(m7-m8+m9-m10+m11-m12+m13-m14)
                + 2.0/105.0*(m19-m20+m21-m22+m23-m24-m25+m26)
                + 5.0/504.0*(m27-m28)
                + 1.0/315.0*(m33-m34+m35-m36+m37-m38+m39-m40+m41-m42+m43-m44+m45-m46+m47-m48)
                + 1.0/630.0*(m57-m58+m59-m60+m61-m62-m63+m64+m65-m66+m67-m68-m69+m70+m71-m72+m73-m74-m75+m76+m77-m78+m79-m80)
                + 1.0/5040.0*(m81-m82+m83-m84+m85-m86+m87-m88);

            ny =  4.0/45.0*(m3-m4)
                + 1.0/21.0*(m7-m8-m9+m10+m15-m16+m17-m18)
                + 2.0/105.0*(m19-m20-m21+m22+m23-m24+m25-m26)
                + 5.0/504.0*(m29-m30)
                + 1.0/315.0*(m33-m34+m35-m36-m37+m38-m39+m40+m49-m50+m51-m52+m53-m54+m55-m56)
                + 1.0/630.0*(m57-m58+m59-m60+m61-m62+m63-m64-m65+m66+m67-m68+m69-m70-m71+m72+m73-m74+m75-m76-m77+m78+m79-m80)
                + 1.0/5040.0*(m81-m82-m83+m84+m89-m90+m91-m92);

            nz =  4.0/45.0*(m5-m6)
                + 1.0/21.0*(m11-m12-m13+m14+m15-m16-m17+m18)
                + 2.0/105.0*(m19-m20+m21-m22-m23+m24+m25-m26)
                + 5.0/504.0*(m31-m32)
                + 1.0/315.0*(m41-m42+m43-m44-m45+m46-m47+m48+m49-m50+m51-m52-m53+m54-m55+m56)
                + 1.0/630.0*(m57-m58+m59-m60+m61-m62+m63-m64+m65-m66-m67+m68+m69-m70+m71-m72-m73+m74+m75-m76+m77-m78-m79+m80)
                + 1.0/5040.0*(m85-m86-m87+m88+m89-m90-m91+m92);

    		//............Compute the Chemical Potential...............................
            //chem = 2.0*3.0/18.0*(m1+m2+m3+m4+m5+m6-6*phi+0.5*(m7+m8+m9+m10+m11+m12+m13+m14+m15+m16+m17+m18-12*phi));//intermediate var, i.e. the laplacian
            //chem = 4.0*beta*phi*(phi+1.0)*(phi-1.0)-kappa*chem;
            chem = 2.0*3.0/18.0*(m1+m2+m3+m4+m5+m6-6*phi_temp+0.5*(m7+m8+m9+m10+m11+m12+m13+m14+m15+m16+m17+m18-12*phi_temp));//intermediate var, i.e. the laplacian
            chem = 4.0*beta*phi_temp*(phi_temp+1.0)*(phi_temp-1.0)-kappa*chem;
    		
            //de-noise color gradient and mixed gradient
            C = sqrt(nx*nx+ny*ny+nz*nz);
            if (C<1.0e-12) nx=ny=nz=0.0;
            //maybe you can also de-noise chemical potential ? within the bulk phase chem should be zero
            if (fabs(chem)<1.0e-12) chem=0.0;

    		// q=0
    		m0 = dist[n];
    		// q=1
    		nr1 = neighborList[n]; // neighbor 2 ( > 10Np => odd part of dist)
    		m1 = dist[nr1]; // reading the f1 data into register fq

    		nr2 = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
    		m2 = dist[nr2];  // reading the f2 data into register fq

    		// q=3
    		nr3 = neighborList[n+2*Np]; // neighbor 4
    		m3 = dist[nr3];

    		// q = 4
    		nr4 = neighborList[n+3*Np]; // neighbor 3
    		m4 = dist[nr4];

    		// q=5
    		nr5 = neighborList[n+4*Np];
    		m5 = dist[nr5];

    		// q = 6
    		nr6 = neighborList[n+5*Np];
    		m6 = dist[nr6];
    		
    		// q=7
    		nr7 = neighborList[n+6*Np];
    		m7 = dist[nr7];

    		// q = 8
    		nr8 = neighborList[n+7*Np];
    		m8 = dist[nr8];

    		// q=9
    		nr9 = neighborList[n+8*Np];
    		m9 = dist[nr9];

    		// q = 10
    		nr10 = neighborList[n+9*Np];
    		m10 = dist[nr10];

    		// q=11
    		nr11 = neighborList[n+10*Np];
    		m11 = dist[nr11];

    		// q=12
    		nr12 = neighborList[n+11*Np];
    		m12 = dist[nr12];

    		// q=13
    		nr13 = neighborList[n+12*Np];
    		m13 = dist[nr13];

    		// q=14
    		nr14 = neighborList[n+13*Np];
    		m14 = dist[nr14];

    		// q=15
    		nr15 = neighborList[n+14*Np];
    		m15 = dist[nr15];

    		// q=16
    		nr16 = neighborList[n+15*Np];
    		m16 = dist[nr16];

    		// q=17
    		nr17 = neighborList[n+16*Np];
    		m17 = dist[nr17];

    		// q=18
    		nr18 = neighborList[n+17*Np];
    		m18 = dist[nr18];

            //compute fluid velocity
            ux = 3.0/rho0*(m1-m2+m7-m8+m9-m10+m11-m12+m13-m14+0.5*(chem*nx+Fx)/3.0);
            uy = 3.0/rho0*(m3-m4+m7-m8-m9+m10+m15-m16+m17-m18+0.5*(chem*ny+Fy)/3.0);
            uz = 3.0/rho0*(m5-m6+m11-m12-m13+m14+m15-m16-m17+m18+0.5*(chem*nz+Fz)/3.0);
            //compute pressure
            p = (m0+m2+m1+m4+m3+m6+m5+m8+m7+m10+m9+m12+m11+m14+m13+m16+m15+m18+m17)
                      +0.5*(rhoA-rhoB)/2.0/3.0*(ux*nx+uy*ny+uz*nz);

            //compute equilibrium distributions
                feq0 = 0.3333333333333333*p - 0.25*(Fx*ux + Fy*uy + Fz*uz)*(-0.6666666666666666 + ux*ux + uy*uy + uz*uz) - 
         0.16666666666666666*rho0*(ux*ux + uy*uy + uz*uz) - 0.5*(-(nx*ux) - ny*uy - nz*uz)*
          (-0.08333333333333333*(rhoA - rhoB)*(ux*ux + uy*uy + uz*uz) + chem*(0.3333333333333333 - 0.5*(ux*ux + uy*uy + uz*uz)));
                feq1 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-ux*ux + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*uz)) - 
         0.125*(Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 
           0.3333333333333333*(-2.*ux + ux*ux + uy*uy + uz*uz)) - 0.0625*(nx - nx*ux - ny*uy - nz*uz)*
          (2*chem*ux*ux - 0.3333333333333333*((-rhoA + rhoB)*ux*ux + 2*chem*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
           0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*uz)));
                feq2 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-ux*ux + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*uz)) - 
         0.125*(Fx + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 
           0.3333333333333333*(2.*ux + ux*ux + uy*uy + uz*uz)) - 0.0625*(nx + nx*ux + ny*uy + nz*uz)*
          (-2.*chem*ux*ux + 0.1111111111111111*(-4.*chem + rhoB*(-2.*ux - 1.*ux*ux - 1.*uy*uy - 1.*uz*uz) + 
             rhoA*(2.*ux + ux*ux + uy*uy + uz*uz)) + 0.3333333333333333*((-1.*rhoA + rhoB)*ux*ux + 
             chem*(4.*ux + 2.*ux*ux + 2.*uy*uy + 2.*uz*uz)));
                feq3 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uy*uy + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*uz)) - 
         0.125*(Fx*ux + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy + 
           0.3333333333333333*(ux*ux - 2.*uy + uy*uy + uz*uz)) - 0.0625*(ny - nx*ux - ny*uy - nz*uz)*
          (2*chem*uy*uy - 0.3333333333333333*((-rhoA + rhoB)*uy*uy + 2*chem*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
           0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*uz)));
                feq4 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uy*uy + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*uz)) - 
         0.125*(Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy + 
           0.3333333333333333*(ux*ux + 2.*uy + uy*uy + uz*uz)) - 0.0625*(ny + nx*ux + ny*uy + nz*uz)*
          (-2.*chem*uy*uy + 0.1111111111111111*(-4.*chem + rhoB*(-1.*ux*ux - 2.*uy - 1.*uy*uy - 1.*uz*uz) + 
             rhoA*(ux*ux + 2.*uy + uy*uy + uz*uz)) + 0.3333333333333333*((-1.*rhoA + rhoB)*uy*uy + 
             chem*(2.*ux*ux + 4.*uy + 2.*uy*uy + 2.*uz*uz))); 
                feq5 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uz*uz + 0.3333333333333333*(ux*ux + uy*uy + (-2 + uz)*uz)) - 
         0.125*(Fx*ux + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*uz*uz + 
           0.3333333333333333*(ux*ux + uy*uy + (-2. + uz)*uz)) - 0.0625*(nx*ux + ny*uy + nz*(-1. + uz))*
          (-2.*chem*uz*uz + 0.1111111111111111*(-4.*chem + rhoB*(-1.*ux*ux - 1.*uy*uy + (2. - 1.*uz)*uz) + 
             rhoA*(ux*ux + uy*uy + (-2. + uz)*uz)) + 0.3333333333333333*((-1.*rhoA + rhoB)*uz*uz + 
             chem*(2.*ux*ux + 2.*uy*uy + uz*(-4. + 2.*uz)))); 
                feq6 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uz*uz + 0.3333333333333333*(ux*ux + uy*uy + uz*(2 + uz))) - 
         0.125*(Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*uz*uz + 
           0.3333333333333333*(ux*ux + uy*uy + uz*(2. + uz))) - 0.0625*(nz + nx*ux + ny*uy + nz*uz)*
          (-2.*chem*uz*uz + 0.1111111111111111*(-4.*chem + rhoB*(-1.*ux*ux - 1.*uy*uy + (-2. - 1.*uz)*uz) + 
             rhoA*(ux*ux + uy*uy + uz*(2. + uz))) + 0.3333333333333333*((-1.*rhoA + rhoB)*uz*uz + 
             chem*(2.*ux*ux + 2.*uy*uy + uz*(4. + 2.*uz)))); 
                feq7 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux + uy)*(ux + uy) + 0.3333333333333333*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) - 
         0.0625*(Fx*(-1. + ux) + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uy - 1.*uy*uy + 
           0.3333333333333333*(-2.*ux + ux*ux - 2.*uy + uy*uy + uz*uz)) - 0.03125*(nx + ny - nx*ux - ny*uy - nz*uz)*
          (2*chem*(ux + uy)*(ux + uy) + 0.3333333333333333*((rhoA - rhoB)*(ux + uy)*(ux + uy) - 2*chem*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
           0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz))); 
                feq8 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux + uy)*(ux + uy) + 0.3333333333333333*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) - 
         0.0625*(Fx + Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uy - 1.*uy*uy + 
           0.3333333333333333*(2.*ux + ux*ux + 2.*uy + uy*uy + uz*uz)) - 0.03125*(-(nx*(1 + ux)) - ny*(1 + uy) - nz*uz)*
          (2*chem*(ux + uy)*(ux + uy) - 0.3333333333333333*(-((rhoA - rhoB)*(ux + uy)*(ux + uy)) + 
             2*chem*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz))); 
                feq9 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux - uy)*(ux - uy) + 0.3333333333333333*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) - 
         0.0625*(Fy + Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uy - 1.*uy*uy + 
           0.3333333333333333*(-2.*ux + ux*ux + 2.*uy + uy*uy + uz*uz)) - 0.03125*(nx - nx*ux - ny*(1 + uy) - nz*uz)*
          (2*chem*(ux - uy)*(ux - uy) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uy)*(ux - uy)) + 
             2*chem*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz))); 
                feq10 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux - uy)*(ux - uy) + 0.3333333333333333*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) - 
         0.0625*(Fx*(1 + ux) + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uy - 1.*uy*uy + 
           0.3333333333333333*(2.*ux + ux*ux - 2.*uy + uy*uy + uz*uz)) - 0.03125*(ny - nx*(1 + ux) - ny*uy - nz*uz)*
          (2*chem*(ux - uy)*(ux - uy) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uy)*(ux - uy)) + 
             2*chem*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz))); 
                feq11 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux + uz)*(ux + uz) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) - 
         0.0625*(Fx*(-1. + ux) + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uz - 1.*uz*uz + 
           0.3333333333333333*(-2.*ux + ux*ux + uy*uy + (-2. + uz)*uz)) - 0.03125*(nx + nz - nx*ux - ny*uy - nz*uz)*
          (2*chem*(ux + uz)*(ux + uz) + 0.3333333333333333*((rhoA - rhoB)*(ux + uz)*(ux + uz) - 2*chem*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
           0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz))); 
                feq12 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux + uz)*(ux + uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) - 
         0.0625*(Fx + Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uz - 1.*uz*uz + 
           0.3333333333333333*(2.*ux + ux*ux + uy*uy + uz*(2. + uz))) - 0.03125*(-(nx*(1 + ux)) - ny*uy - nz*(1 + uz))*
          (2*chem*(ux + uz)*(ux + uz) - 0.3333333333333333*(-((rhoA - rhoB)*(ux + uz)*(ux + uz)) + 
             2*chem*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))); 
                feq13 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux - uz)*(ux - uz) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) - 
         0.0625*(Fz + Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uz - 1.*uz*uz + 
           0.3333333333333333*(-2.*ux + ux*ux + uy*uy + uz*(2. + uz))) - 0.03125*(nx - nx*ux - ny*uy - nz*(1 + uz))*
          (2*chem*(ux - uz)*(ux - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uz)*(ux - uz)) + 
             2*chem*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*(2 + uz)))); 
                feq14 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux - uz)*(ux - uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) - 
         0.0625*(Fx*(1 + ux) + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uz - 1.*uz*uz + 
           0.3333333333333333*(2.*ux + ux*ux + uy*uy + (-2. + uz)*uz)) - 0.03125*(nz - nx*(1 + ux) - ny*uy - nz*uz)*
          (2*chem*(ux - uz)*(ux - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uz)*(ux - uz)) + 
             2*chem*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz))); 
                feq15 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) - 
         0.0625*(Fx*ux + Fy*(-1. + uy) + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*uy*uy - 2.*uy*uz - 1.*uz*uz + 
           0.3333333333333333*(ux*ux - 2.*uy + uy*uy + (-2. + uz)*uz)) - 0.03125*(ny + nz - nx*ux - ny*uy - nz*uz)*
          (2*chem*(uy + uz)*(uy + uz) + 0.3333333333333333*((rhoA - rhoB)*(uy + uz)*(uy + uz) - 2*chem*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) + 
           0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))); 
                feq16 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) - 
         0.0625*(Fy + Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy - 2.*uy*uz - 1.*uz*uz + 
           0.3333333333333333*(ux*ux + 2.*uy + uy*uy + uz*(2. + uz))) - 0.03125*(-(nx*ux) - ny*(1 + uy) - nz*(1 + uz))*
          (2*chem*(uy + uz)*(uy + uz) - 0.3333333333333333*(-((rhoA - rhoB)*(uy + uz)*(uy + uz)) + 
             2*chem*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))); 
                feq17 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) - 
         0.0625*(Fz + Fx*ux + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy + 2.*uy*uz - 1.*uz*uz + 
           0.3333333333333333*(ux*ux - 2.*uy + uy*uy + uz*(2. + uz))) - 0.03125*(ny - nx*ux - ny*uy - nz*(1 + uz))*
          (2*chem*(uy - uz)*(uy - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(uy - uz)*(uy - uz)) + 
             2*chem*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))); 
                feq18 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) - 
         0.0625*(Fx*ux + Fy*(1 + uy) + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*uy*uy + 2.*uy*uz - 1.*uz*uz + 
           0.3333333333333333*(ux*ux + 2.*uy + uy*uy + (-2. + uz)*uz)) - 0.03125*(nz - nx*ux - ny*(1 + uy) - nz*uz)*
          (2*chem*(uy - uz)*(uy - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(uy - uz)*(uy - uz)) + 
             2*chem*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz))); 

            //------------------------------------------------- BCK collison ------------------------------------------------------------//
    		// q=0
    		dist[n] = m0 - (m0-feq0)/tau + 0.25*(2*(Fx*ux + Fy*uy + Fz*uz)*(-0.6666666666666666 + ux*ux + uy*uy + uz*uz) + 
         (nx*ux + ny*uy + nz*uz)*(2*chem*(ux*ux + uy*uy + uz*uz) + 
            0.3333333333333333*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + uz*uz)))); 

    		// q = 1
    		dist[nr2] = m1 - (m1-feq1)/tau + 0.125*(2*(Fx*(-1 + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - ux*ux + 
            0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
         (nx*(-1 + ux) + ny*uy + nz*uz)*(-2*chem*(ux*ux) + 
            0.3333333333333333*((-rhoA + rhoB)*(ux*ux) + 2*chem*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*uz))));

    		// q=2
    		dist[nr1] = m2 - (m2-feq2)/tau + 0.125*(2*(Fx + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - ux*ux + 
            0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*uz)) + 
         (nx + nx*ux + ny*uy + nz*uz)*(-2*chem*(ux*ux) + 
            0.3333333333333333*((-rhoA + rhoB)*(ux*ux) + 2*chem*(2*ux + ux*ux + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*uz))));

    		// q = 3
    		dist[nr4] = m3 - (m3-feq3)/tau + 0.125*(2*(Fx*ux + Fy*(-1 + uy) + Fz*uz)*(-0.2222222222222222 - uy*uy + 
            0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
         (nx*ux + ny*(-1 + uy) + nz*uz)*(-2*chem*(uy*uy) + 
            0.3333333333333333*((-rhoA + rhoB)*(uy*uy) + 2*chem*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*uz))));

    		// q = 4
    		dist[nr3] = m4 - (m4-feq4)/tau + 0.125*(2*(Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - uy*uy + 
            0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*uz)) + 
         (ny + nx*ux + ny*uy + nz*uz)*(-2*chem*(uy*uy) + 
            0.3333333333333333*((-rhoA + rhoB)*(uy*uy) + 2*chem*(ux*ux + 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*uz))));

    		// q = 5
    		dist[nr6] = m5 - (m5-feq5)/tau + 0.125*(2*(Fx*ux + Fy*uy + Fz*(-1 + uz))*(-0.2222222222222222 - uz*uz + 
            0.3333333333333333*(ux*ux + uy*uy + (-2 + uz)*uz)) + 
         (nx*ux + ny*uy + nz*(-1 + uz))*(-2*chem*(uz*uz) + 
            0.3333333333333333*((-rhoA + rhoB)*(uz*uz) + 2*chem*(ux*ux + uy*uy + (-2 + uz)*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + (-2 + uz)*uz))));

    		// q = 6
    		dist[nr5] = m6 - (m6-feq6)/tau + 0.125*(2*(Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - uz*uz + 
            0.3333333333333333*(ux*ux + uy*uy + uz*(2 + uz))) + 
         (nz + nx*ux + ny*uy + nz*uz)*(-2*chem*(uz*uz) + 
            0.3333333333333333*((-rhoA + rhoB)*(uz*uz) + 2*chem*(ux*ux + uy*uy + uz*(2 + uz))) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + uz*(2 + uz)))));

    		// q = 7
    		dist[nr8] = m7 - (m7-feq7)/tau + 0.0625*(-2*(Fx*(-1 + ux) + Fy*(-1 + uy) + Fz*uz)*
          (0.2222222222222222 + (ux + uy)*(ux + uy) - 
            0.3333333333333333*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
         (nx*(-1 + ux) + ny*(-1 + uy) + nz*uz)*
          (-2*chem*((ux + uy)*(ux + uy)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux + uy)*(ux + uy))) + 2*chem*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz))));

    		// q = 8
    		dist[nr7] = m8 - (m8-feq8)/tau + 0.0625*(2*(Fx + Fy + Fx*ux + Fy*uy + Fz*uz)*
          (-0.2222222222222222 - (ux + uy)*(ux + uy) + 
            0.3333333333333333*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
         (nx + ny + nx*ux + ny*uy + nz*uz)*
          (-2*chem*((ux + uy)*(ux + uy)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux + uy)*(ux + uy))) + 2*chem*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz))));

    		// q = 9
    		dist[nr10] = m9 - (m9-feq9)/tau + 0.0625*(2*(Fy + Fx*(-1 + ux) + Fy*uy + Fz*uz)*
          (-0.2222222222222222 - (ux - uy)*(ux - uy) + 
            0.3333333333333333*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
         (ny + nx*(-1 + ux) + ny*uy + nz*uz)*
          (-2*chem*((ux - uy)*(ux - uy)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux - uy)*(ux - uy))) + 2*chem*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz))));

    		// q = 10
    		dist[nr9] = m10 - (m10-feq10)/tau + 0.0625*(2*(Fx*(1 + ux) + Fy*(-1 + uy) + Fz*uz)*
          (-0.2222222222222222 - (ux - uy)*(ux - uy) + 
            0.3333333333333333*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
         (nx*(1 + ux) + ny*(-1 + uy) + nz*uz)*
          (-2*chem*((ux - uy)*(ux - uy)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux - uy)*(ux - uy))) + 2*chem*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz))));

    		// q = 11
    		dist[nr12] = m11 - (m11-feq11)/tau + 0.0625*(-2*(Fx*(-1 + ux) + Fy*uy + Fz*(-1 + uz))*
          (0.2222222222222222 + (ux + uz)*(ux + uz) - 
            0.3333333333333333*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
         (nx*(-1 + ux) + ny*uy + nz*(-1 + uz))*
          (-2*chem*((ux + uz)*(ux + uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux + uz)*(ux + uz))) + 2*chem*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz))));

    		// q = 12
    		dist[nr11] = m12 - (m12-feq12)/tau + 0.0625*(2*(Fx + Fz + Fx*ux + Fy*uy + Fz*uz)*
          (-0.2222222222222222 - (ux + uz)*(ux + uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))
           + (nx + nz + nx*ux + ny*uy + nz*uz)*
          (-2*chem*((ux + uz)*(ux + uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux + uz)*(ux + uz))) + 2*chem*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))));

    		// q = 13
    		dist[nr14] = m13 - (m13-feq13)/tau + 0.0625*(2*(Fz + Fx*(-1 + ux) + Fy*uy + Fz*uz)*
          (-0.2222222222222222 - (ux - uz)*(ux - uz) + 
            0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
         (nz + nx*(-1 + ux) + ny*uy + nz*uz)*
          (-2*chem*((ux - uz)*(ux - uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux - uz)*(ux - uz))) + 2*chem*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*(2 + uz)))));

    		// q= 14
    		dist[nr13] = m14 - (m14-feq14)/tau + 0.0625*(2*(Fx*(1 + ux) + Fy*uy + Fz*(-1 + uz))*
          (-0.2222222222222222 - (ux - uz)*(ux - uz) + 
            0.3333333333333333*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
         (nx*(1 + ux) + ny*uy + nz*(-1 + uz))*
          (-2*chem*((ux - uz)*(ux - uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux - uz)*(ux - uz))) + 2*chem*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz))));

    		// q = 15
    		dist[nr16] = m15 - (m15-feq15)/tau + 0.0625*(-2*(Fx*ux + Fy*(-1 + uy) + Fz*(-1 + uz))*
          (0.2222222222222222 + (uy + uz)*(uy + uz) - 0.3333333333333333*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))
           + (nx*ux + ny*(-1 + uy) + nz*(-1 + uz))*
          (-2*chem*((uy + uz)*(uy + uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((uy + uz)*(uy + uz))) + 2*chem*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))));

    		// q = 16
    		dist[nr15] = m16 - (m16-feq16)/tau + 0.0625*(2*(Fy + Fz + Fx*ux + Fy*uy + Fz*uz)*
          (-0.2222222222222222 - (uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))
           + (ny + nz + nx*ux + ny*uy + nz*uz)*
          (-2*chem*((uy + uz)*(uy + uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((uy + uz)*(uy + uz))) + 2*chem*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))));

    		// q = 17
    		dist[nr18] = m17 - (m17-feq17)/tau + 0.0625*(2*(Fz + Fx*ux + Fy*(-1 + uy) + Fz*uz)*
          (-0.2222222222222222 - (uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))
           + (nz + nx*ux + ny*(-1 + uy) + nz*uz)*
          (-2*chem*((uy - uz)*(uy - uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((uy - uz)*(uy - uz))) + 2*chem*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))));

    		// q = 18
    		dist[nr17] = m18 - (m18-feq18)/tau + 0.0625*(2*(Fx*ux + Fy*(1 + uy) + Fz*(-1 + uz))*
          (-0.2222222222222222 - (uy - uz)*(uy - uz) + 
            0.3333333333333333*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 
         (nx*ux + ny*(1 + uy) + nz*(-1 + uz))*
          (-2*chem*((uy - uz)*(uy - uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((uy - uz)*(uy - uz))) + 2*chem*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz))));
            //----------------------------------------------------------------------------------------------------------------------------------------//

            // ----------------------------- compute phase field evolution ----------------------------------------
            //Normalize the Color Gradient
            C = sqrt(nx*nx+ny*ny+nz*nz);
            double ColorMag = C;
            if (C==0.0) ColorMag=1.0;
            nx = nx/ColorMag;
            ny = ny/ColorMag;
            nz = nz/ColorMag;		
            //compute surface tension-related parameter
            //theta = 4.5*M*2.0*(1-phi*phi)/W;
            theta = 4.5*M*2.0*(1-phi_temp*phi_temp)/W;

            //load distributions of phase field
            //q=0
            h0 = hq[n];
            //q=1
            h1 = hq[nr1]; 

            //q=2
            h2 = hq[nr2];  

            //q=3
            h3 = hq[nr3];

            //q=4
            h4 = hq[nr4];

            //q=5
            h5 = hq[nr5];

            //q=6
            h6 = hq[nr6];

            //-------------------------------- BGK collison for phase field ---------------------------------//
            // q = 0
            hq[n] = h0 - (h0 - 0.3333333333333333*phi)/tauM;

            // q = 1
            hq[nr2] = h1 - (h1 - 0.1111111111111111*nx*theta - phi*(0.1111111111111111 + 0.5*ux))/tauM;

            // q = 2
            hq[nr1] = h2 - (h2 + 0.1111111111111111*nx*theta - phi*(0.1111111111111111 - 0.5*ux))/tauM;

            // q = 3
            hq[nr4] = h3 - (h3 - 0.1111111111111111*ny*theta - phi*(0.1111111111111111 + 0.5*uy))/tauM;

            // q = 4
            hq[nr3] = h4 - (h4 + 0.1111111111111111*ny*theta - phi*(0.1111111111111111 - 0.5*uy))/tauM;

            // q = 5
            hq[nr6] = h5 - (h5 - 0.1111111111111111*nz*theta - phi*(0.1111111111111111 + 0.5*uz))/tauM;

            // q = 6
            hq[nr5] = h6 - (h6 + 0.1111111111111111*nz*theta - phi*(0.1111111111111111 - 0.5*uz))/tauM;
            //........................................................................

            //Update velocity on device
    		Vel[0*Np+n] = ux;
    		Vel[1*Np+n] = uy;
    		Vel[2*Np+n] = uz;
            //Update pressure on device
            Pressure[n] = p;
            //Update chemical potential on device
            mu_phi[n] = chem;
            //Update color gradient on device
    		ColorGrad[0*Np+n] = nx;
    		ColorGrad[1*Np+n] = ny;
    		ColorGrad[2*Np+n] = nz;
    	}
    }
}

__global__ void dvc_ScaLBL_D3Q19_AAeven_FreeLeeModel_Combined_HigherOrder(int *Map, double *dist, double *hq, double *Den,	double *Phi, double *mu_phi, double *Vel, double *Pressure, double *ColorGrad,
        double rhoA, double rhoB, double tauA, double tauB, double tauM, double kappa, double beta, double W, double Fx, double Fy, double Fz, 
        int strideY, int strideZ, int start, int finish, int Np){

	int n,nn,ijk;
    double ux,uy,uz;//fluid velocity 
    double p;//pressure
    double chem;//chemical potential
    double phi; //phase field
    double rho0;//fluid density
	// register variables for neighbors of phase field
    double m0;
    // w(|c|^2) = 1 
	double m1,m2,m3,m4,m5,m6;
    // w(|c|^2) = 2
	double m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
    // w(|c|^2) = 3
    double m19,m20,m21,m22,m23,m24,m25,m26;
    // w(|c|^2) = 4
    double m27,m28,m29,m30,m31,m32;
    // w(|c|^2) = 5
    double m33,m34,m35,m36,m37,m38,m39,m40,m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,m51,m52,m53,m54,m55,m56;
    // w(|c|^2) = 6
    double m57,m58,m59,m60,m61,m62,m63,m64,m65,m66,m67,m68,m69,m70,m71,m72,m73,m74,m75,m76,m77,m78,m79,m80;
    // w(|c|^2) = 8
    double m81,m82,m83,m84,m85,m86,m87,m88,m89,m90,m91,m92;
	//double mm1,mm2,mm4,mm6,mm8,mm9,mm10,mm11,mm12,mm13,mm14,mm15,mm16,mm17,mm18;
	//double mm3,mm5,mm7;
    double feq0,feq1,feq2,feq3,feq4,feq5,feq6,feq7,feq8,feq9,feq10,feq11,feq12,feq13,feq14,feq15,feq16,feq17,feq18;
    double nx,ny,nz;//normal color gradient
    //double mgx,mgy,mgz;//mixed gradient reaching secondary neighbor

    double h0,h1,h2,h3,h4,h5,h6;//distributions for LB phase field
	double tau;//position dependent LB relaxation time for fluid
    double C,theta;
    double M = 2.0/9.0*(tauM-0.5);//diffusivity (or mobility) for the phase field D3Q7
    double phi_temp;

    int S = Np/NBLOCKS/NTHREADS + 1;
    for (int s=0; s<S; s++){
    	//........Get 1-D index for this thread....................
    	n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

    	if ( n<finish ){
    		rho0 = Den[n];//load density

    		// Get the 1D index based on regular data layout
    		ijk = Map[n];
            phi = Phi[ijk];// load phase field
            phi_temp = phi;
            if (phi>1.f) phi_temp=1.0;
            if (phi<-1.f) phi_temp=-1.0;

    		// local relaxation time
    		tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);

    		//					COMPUTE THE COLOR GRADIENT
    		//........................................................................
    		//.................Read Phase Indicator Values............................
    		//........................................................................
            //-------------------
            //--- w(|c|^2) = 1--- 
            //-------------------
    		nn = ijk-1;							// neighbor index (get convention)
    		m2 = Phi[nn];						// get neighbor for phi - 1
    		//........................................................................
    		nn = ijk+1;							// neighbor index (get convention)
    		m1 = Phi[nn];						// get neighbor for phi - 2
    		//........................................................................
    		nn = ijk-strideY;							// neighbor index (get convention)
    		m4 = Phi[nn];					// get neighbor for phi - 3
    		//........................................................................
    		nn = ijk+strideY;							// neighbor index (get convention)
    		m3 = Phi[nn];					// get neighbor for phi - 4
    		//........................................................................
    		nn = ijk-strideZ;						// neighbor index (get convention)
    		m6 = Phi[nn];					// get neighbor for phi - 5
    		//........................................................................
    		nn = ijk+strideZ;						// neighbor index (get convention)
    		m5 = Phi[nn];					// get neighbor for phi - 6
    		//........................................................................

            //-------------------
            //--- w(|c|^2) = 2--- 
            //-------------------
    		nn = ijk-strideY-1;						// neighbor index (get convention)
    		m8 = Phi[nn];					// get neighbor for phi - 7
    		//........................................................................
    		nn = ijk+strideY+1;						// neighbor index (get convention)
    		m7 = Phi[nn];					// get neighbor for phi - 8
    		//........................................................................
    		nn = ijk+strideY-1;						// neighbor index (get convention)
    		m10 = Phi[nn];					// get neighbor for phi - 9
    		//........................................................................
    		nn = ijk-strideY+1;						// neighbor index (get convention)
    		m9 = Phi[nn];					// get neighbor for phi - 10
    		//........................................................................
    		nn = ijk-strideZ-1;						// neighbor index (get convention)
    		m12 = Phi[nn];					// get neighbor for phi - 11
    		//........................................................................
    		nn = ijk+strideZ+1;						// neighbor index (get convention)
    		m11 = Phi[nn];					// get neighbor for phi - 12
    		//........................................................................
    		nn = ijk+strideZ-1;						// neighbor index (get convention)
    		m14 = Phi[nn];					// get neighbor for phi - 13
    		//........................................................................
    		nn = ijk-strideZ+1;						// neighbor index (get convention)
    		m13 = Phi[nn];					// get neighbor for phi - 14
    		//........................................................................
    		nn = ijk-strideZ-strideY;					// neighbor index (get convention)
    		m16 = Phi[nn];					// get neighbor for phi - 15
    		//........................................................................
    		nn = ijk+strideZ+strideY;					// neighbor index (get convention)
    		m15 = Phi[nn];					// get neighbor for phi - 16
    		//........................................................................
    		nn = ijk+strideZ-strideY;					// neighbor index (get convention)
    		m18 = Phi[nn];					// get neighbor for phi - 17
    		//........................................................................
    		nn = ijk-strideZ+strideY;					// neighbor index (get convention)
    		m17 = Phi[nn];					// get neighbor for phi - 18
    		//........................................................................

            //-------------------
            //--- w(|c|^2) = 3--- 
            //-------------------
            // c19 = (1,1,1) = (cx,cy,cz)
    		nn = ijk+strideZ+strideY+1;	
    		m19 = Phi[nn];			
    		//........................................................................
            // c20 = (-1,-1,-1)
    		nn = ijk-strideZ-strideY-1;	
    		m20 = Phi[nn];			
    		//........................................................................
            // c21 = (1,-1,1)
    		nn = ijk+strideZ-strideY+1;	
    		m21 = Phi[nn];			
    		//........................................................................
            // c22 = (-1,1,-1)
    		nn = ijk-strideZ+strideY-1;	
    		m22 = Phi[nn];			
    		//........................................................................
            // c23 = (1,1,-1)
    		nn = ijk-strideZ+strideY+1;	
    		m23 = Phi[nn];			
    		//........................................................................
            // c24 = (-1,-1,1)
    		nn = ijk+strideZ-strideY-1;	
    		m24 = Phi[nn];			
    		//........................................................................
            // c25 = (-1,1,1)
    		nn = ijk+strideZ+strideY-1;	
    		m25 = Phi[nn];			
    		//........................................................................
            // c26 = (1,-1,-1)
    		nn = ijk-strideZ-strideY+1;	
    		m26 = Phi[nn];			
    		//........................................................................
            
            //-------------------
            //--- w(|c|^2) = 4--- 
            //-------------------
            // c27 = (2,0,0) = (cx,cy,cz)
    		nn = ijk+2;		
    		m27 = Phi[nn];			
    		//........................................................................
            // c28 = (-2,0,0) = (cx,cy,cz)
    		nn = ijk-2;		
    		m28 = Phi[nn];			
    		//........................................................................
            // c29 = (0,2,0) = (cx,cy,cz)
    		nn = ijk+2*strideY;		
    		m29 = Phi[nn];			
    		//........................................................................
            // c30 = (0,-2,0) = (cx,cy,cz)
    		nn = ijk-2*strideY;		
    		m30 = Phi[nn];			
    		//........................................................................
            // c31 = (0,0,2) = (cx,cy,cz)
    		nn = ijk+2*strideZ;		
    		m31 = Phi[nn];			
    		//........................................................................
            // c32 = (0,0,-2) = (cx,cy,cz)
    		nn = ijk-2*strideZ;		
    		m32 = Phi[nn];			
    		//........................................................................

            //-------------------
            //--- w(|c|^2) = 5--- 
            //-------------------
            // c33 = (1,2,0) = (cx,cy,cz)
    		nn = ijk+2*strideY+1;				
    		m33 = Phi[nn];					
    		//........................................................................
            // c34 = (-1,-2,0) = (cx,cy,cz)
    		nn = ijk-2*strideY-1;				
    		m34 = Phi[nn];					
    		//........................................................................
            // c35 = (2,1,0) = (cx,cy,cz)
    		nn = ijk+strideY+2;				
    		m35 = Phi[nn];					
    		//........................................................................
            // c36 = (-2,-1,0) = (cx,cy,cz)
    		nn = ijk-strideY-2;				
    		m36 = Phi[nn];					
    		//........................................................................
            // c37 = (1,-2,0) = (cx,cy,cz)
    		nn = ijk-2*strideY+1;				
    		m37 = Phi[nn];					
    		//........................................................................
            // c38 = (-1,2,0) = (cx,cy,cz)
    		nn = ijk+2*strideY-1;				
    		m38 = Phi[nn];					
    		//........................................................................
            // c39 = (2,-1,0) = (cx,cy,cz)
    		nn = ijk-strideY+2;				
    		m39 = Phi[nn];					
    		//........................................................................
            // c40 = (-2,1,0) = (cx,cy,cz)
    		nn = ijk+strideY-2;				
    		m40 = Phi[nn];					
    		//........................................................................
            // c41 = (1,0,2) = (cx,cy,cz)
    		nn = ijk+2*strideZ+1;				
    		m41 = Phi[nn];					
    		//........................................................................
            // c42 = (-1,0,-2) = (cx,cy,cz)
    		nn = ijk-2*strideZ-1;				
    		m42 = Phi[nn];					
    		//........................................................................
            // c43 = (2,0,1) = (cx,cy,cz)
    		nn = ijk+strideZ+2;				
    		m43 = Phi[nn];					
    		//........................................................................
            // c44 = (-2,0,-1) = (cx,cy,cz)
    		nn = ijk-strideZ-2;				
    		m44 = Phi[nn];					
    		//........................................................................
            // c45 = (1,0,-2) = (cx,cy,cz)
    		nn = ijk-2*strideZ+1;				
    		m45 = Phi[nn];					
    		//........................................................................
            // c46 = (-1,0,2) = (cx,cy,cz)
    		nn = ijk+2*strideZ-1;				
    		m46 = Phi[nn];					
    		//........................................................................
            // c47 = (2,0,-1) = (cx,cy,cz)
    		nn = ijk-strideZ+2;				
    		m47 = Phi[nn];					
    		//........................................................................
            // c48 = (-2,0,1) = (cx,cy,cz)
    		nn = ijk+strideZ-2;				
    		m48 = Phi[nn];					
    		//........................................................................
            // c49 = (0,1,2) = (cx,cy,cz)
    		nn = ijk+2*strideZ+strideY;				
    		m49 = Phi[nn];					
    		//........................................................................
            // c50 = (0,-1,-2) = (cx,cy,cz)
    		nn = ijk-2*strideZ-strideY;				
    		m50 = Phi[nn];					
    		//........................................................................
            // c51 = (0,2,1) = (cx,cy,cz)
    		nn = ijk+strideZ+2*strideY;				
    		m51 = Phi[nn];					
    		//........................................................................
            // c52 = (0,-2,-1) = (cx,cy,cz)
    		nn = ijk-strideZ-2*strideY;				
    		m52 = Phi[nn];					
    		//........................................................................
            // c53 = (0,1,-2) = (cx,cy,cz)
    		nn = ijk-2*strideZ+strideY;				
    		m53 = Phi[nn];					
    		//........................................................................
            // c54 = (0,-1,2) = (cx,cy,cz)
    		nn = ijk+2*strideZ-strideY;				
    		m54 = Phi[nn];					
    		//........................................................................
            // c55 = (0,2,-1) = (cx,cy,cz)
    		nn = ijk-strideZ+2*strideY;				
    		m55 = Phi[nn];					
    		//........................................................................
            // c56 = (0,-2,1) = (cx,cy,cz)
    		nn = ijk+strideZ-2*strideY;				
    		m56 = Phi[nn];					
    		//........................................................................

            //-------------------
            //--- w(|c|^2) = 6--- 
            //-------------------
            // c57 = (1,1,2) = (cx,cy,cz)
    		nn = ijk+1+strideY+2*strideZ;				
    		m57 = Phi[nn];					
    		//........................................................................
            // c58 = (-1,-1,-2) = (cx,cy,cz)
    		nn = ijk-1-strideY-2*strideZ;				
    		m58 = Phi[nn];					
    		//........................................................................
            // c59 = (1,2,1) = (cx,cy,cz)
    		nn = ijk+1+2*strideY+strideZ;				
    		m59 = Phi[nn];					
    		//........................................................................
            // c60 = (-1,-2,-1) = (cx,cy,cz)
    		nn = ijk-1-2*strideY-strideZ;				
    		m60 = Phi[nn];					
    		//........................................................................
            // c61 = (2,1,1) = (cx,cy,cz)
    		nn = ijk+2+strideY+strideZ;				
    		m61 = Phi[nn];					
    		//........................................................................
            // c62 = (-2,-1,-1) = (cx,cy,cz)
    		nn = ijk-2-strideY-strideZ;				
    		m62 = Phi[nn];					
    		//........................................................................
            // c63 = (-1,1,2) = (cx,cy,cz)
    		nn = ijk-1+strideY+2*strideZ;				
    		m63 = Phi[nn];					
    		//........................................................................
            // c64 = (1,-1,-2) = (cx,cy,cz)
    		nn = ijk+1-strideY-2*strideZ;				
    		m64 = Phi[nn];					
    		//........................................................................
            // c65 = (1,-1,2) = (cx,cy,cz)
    		nn = ijk+1-strideY+2*strideZ;				
    		m65 = Phi[nn];					
    		//........................................................................
            // c66 = (-1,1,-2) = (cx,cy,cz)
    		nn = ijk-1+strideY-2*strideZ;				
    		m66 = Phi[nn];					
    		//........................................................................
            // c67 = (1,1,-2) = (cx,cy,cz)
    		nn = ijk+1+strideY-2*strideZ;				
    		m67 = Phi[nn];					
    		//........................................................................
            // c68 = (-1,-1,2) = (cx,cy,cz)
    		nn = ijk-1-strideY+2*strideZ;				
    		m68 = Phi[nn];					
    		//........................................................................
            // c69 = (-1,2,1) = (cx,cy,cz)
    		nn = ijk-1+2*strideY+strideZ;				
    		m69 = Phi[nn];					
    		//........................................................................
            // c70 = (1,-2,-1) = (cx,cy,cz)
    		nn = ijk+1-2*strideY-strideZ;				
    		m70 = Phi[nn];					
    		//........................................................................
            // c71 = (1,-2,1) = (cx,cy,cz)
    		nn = ijk+1-2*strideY+strideZ;				
    		m71 = Phi[nn];					
    		//........................................................................
            // c72 = (-1,2,-1) = (cx,cy,cz)
    		nn = ijk-1+2*strideY-strideZ;				
    		m72 = Phi[nn];					
    		//........................................................................
            // c73 = (1,2,-1) = (cx,cy,cz)
    		nn = ijk+1+2*strideY-strideZ;				
    		m73 = Phi[nn];					
    		//........................................................................
            // c74 = (-1,-2,1) = (cx,cy,cz)
    		nn = ijk-1-2*strideY+strideZ;				
    		m74 = Phi[nn];					
    		//........................................................................
            // c75 = (-2,1,1) = (cx,cy,cz)
    		nn = ijk-2+strideY+strideZ;				
    		m75 = Phi[nn];					
    		//........................................................................
            // c76 = (2,-1,-1) = (cx,cy,cz)
    		nn = ijk+2-strideY-strideZ;				
    		m76 = Phi[nn];					
    		//........................................................................
            // c77 = (2,-1,1) = (cx,cy,cz)
    		nn = ijk+2-strideY+strideZ;				
    		m77 = Phi[nn];					
    		//........................................................................
            // c78 = (-2,1,-1) = (cx,cy,cz)
    		nn = ijk-2+strideY-strideZ;				
    		m78 = Phi[nn];					
    		//........................................................................
            // c79 = (2,1,-1) = (cx,cy,cz)
    		nn = ijk+2+strideY-strideZ;				
    		m79 = Phi[nn];					
    		//........................................................................
            // c80 = (-2,-1,1) = (cx,cy,cz)
    		nn = ijk-2-strideY+strideZ;				
    		m80 = Phi[nn];					
    		//........................................................................

            //-------------------
            //--- w(|c|^2) = 8--- 
            //-------------------
            // c81 = (2,2,0) = (cx,cy,cz)
    		nn = ijk+2+2*strideY;				
    		m81 = Phi[nn];					
    		//........................................................................
            // c82 = (-2,-2,0) = (cx,cy,cz)
    		nn = ijk-2-2*strideY;				
    		m82 = Phi[nn];					
    		//........................................................................
            // c83 = (2,-2,0) = (cx,cy,cz)
    		nn = ijk+2-2*strideY;				
    		m83 = Phi[nn];					
    		//........................................................................
            // c84 = (-2,2,0) = (cx,cy,cz)
    		nn = ijk-2+2*strideY;				
    		m84 = Phi[nn];					
    		//........................................................................
            // c85 = (2,0,2) = (cx,cy,cz)
    		nn = ijk+2+2*strideZ;				
    		m85 = Phi[nn];					
    		//........................................................................
            // c86 = (-2,0,-2) = (cx,cy,cz)
    		nn = ijk-2-2*strideZ;				
    		m86 = Phi[nn];					
    		//........................................................................
            // c87 = (2,0,-2) = (cx,cy,cz)
    		nn = ijk+2-2*strideZ;				
    		m87 = Phi[nn];					
    		//........................................................................
            // c88 = (-2,0,2) = (cx,cy,cz)
    		nn = ijk-2+2*strideZ;				
    		m88 = Phi[nn];					
    		//........................................................................
            // c89 = (0,2,2) = (cx,cy,cz)
    		nn = ijk+2*strideY+2*strideZ;				
    		m89 = Phi[nn];					
    		//........................................................................
            // c90 = (0,-2,-2) = (cx,cy,cz)
    		nn = ijk-2*strideY-2*strideZ;				
    		m90 = Phi[nn];					
    		//........................................................................
            // c91 = (0,2,-2) = (cx,cy,cz)
    		nn = ijk+2*strideY-2*strideZ;				
    		m91 = Phi[nn];					
    		//........................................................................
            // c92 = (0,-2,2) = (cx,cy,cz)
    		nn = ijk-2*strideY+2*strideZ;				
    		m92 = Phi[nn];					
    		//........................................................................

    		//............Compute the Color Gradient...................................
            nx =  4.0/45.0*(m1-m2)
                + 1.0/21.0*(m7-m8+m9-m10+m11-m12+m13-m14)
                + 2.0/105.0*(m19-m20+m21-m22+m23-m24-m25+m26)
                + 5.0/504.0*(m27-m28)
                + 1.0/315.0*(m33-m34+m35-m36+m37-m38+m39-m40+m41-m42+m43-m44+m45-m46+m47-m48)
                + 1.0/630.0*(m57-m58+m59-m60+m61-m62-m63+m64+m65-m66+m67-m68-m69+m70+m71-m72+m73-m74-m75+m76+m77-m78+m79-m80)
                + 1.0/5040.0*(m81-m82+m83-m84+m85-m86+m87-m88);

            ny =  4.0/45.0*(m3-m4)
                + 1.0/21.0*(m7-m8-m9+m10+m15-m16+m17-m18)
                + 2.0/105.0*(m19-m20-m21+m22+m23-m24+m25-m26)
                + 5.0/504.0*(m29-m30)
                + 1.0/315.0*(m33-m34+m35-m36-m37+m38-m39+m40+m49-m50+m51-m52+m53-m54+m55-m56)
                + 1.0/630.0*(m57-m58+m59-m60+m61-m62+m63-m64-m65+m66+m67-m68+m69-m70-m71+m72+m73-m74+m75-m76-m77+m78+m79-m80)
                + 1.0/5040.0*(m81-m82-m83+m84+m89-m90+m91-m92);

            nz =  4.0/45.0*(m5-m6)
                + 1.0/21.0*(m11-m12-m13+m14+m15-m16-m17+m18)
                + 2.0/105.0*(m19-m20+m21-m22-m23+m24+m25-m26)
                + 5.0/504.0*(m31-m32)
                + 1.0/315.0*(m41-m42+m43-m44-m45+m46-m47+m48+m49-m50+m51-m52-m53+m54-m55+m56)
                + 1.0/630.0*(m57-m58+m59-m60+m61-m62+m63-m64+m65-m66-m67+m68+m69-m70+m71-m72-m73+m74+m75-m76+m77-m78-m79+m80)
                + 1.0/5040.0*(m85-m86-m87+m88+m89-m90-m91+m92);

    		//............Compute the Chemical Potential...............................
            //chem = 2.0*3.0/18.0*(m1+m2+m3+m4+m5+m6-6*phi+0.5*(m7+m8+m9+m10+m11+m12+m13+m14+m15+m16+m17+m18-12*phi));//intermediate var, i.e. the laplacian
            //chem = 4.0*beta*phi*(phi+1.0)*(phi-1.0)-kappa*chem;
            chem = 2.0*3.0/18.0*(m1+m2+m3+m4+m5+m6-6*phi_temp+0.5*(m7+m8+m9+m10+m11+m12+m13+m14+m15+m16+m17+m18-12*phi_temp));//intermediate var, i.e. the laplacian
            chem = 4.0*beta*phi_temp*(phi_temp+1.0)*(phi_temp-1.0)-kappa*chem;

            //de-noise color gradient and mixed gradient
            C = sqrt(nx*nx+ny*ny+nz*nz);
            if (C<1.0e-12) nx=ny=nz=0.0;
            //maybe you can also de-noise chemical potential ? within the bulk phase chem should be ZERO
            if (fabs(chem)<1.0e-12) chem=0.0;
    		
    		// q=0
    		m0 = dist[n];
    		// q=1
    		m1 = dist[2*Np+n]; 

            // q=2
    		m2 = dist[1*Np+n];  

    		// q=3
    		m3 = dist[4*Np+n];

    		// q = 4
    		m4 = dist[3*Np+n];

    		// q=5
    		m5 = dist[6*Np+n];

    		// q = 6
    		m6 = dist[5*Np+n];
    		
    		// q=7
    		m7 = dist[8*Np+n];

    		// q = 8
    		m8 = dist[7*Np+n];

    		// q=9
    		m9 = dist[10*Np+n];

    		// q = 10
    		m10 = dist[9*Np+n];

    		// q=11
    		m11 = dist[12*Np+n];

    		// q=12
    		m12 = dist[11*Np+n];

    		// q=13
    		m13 = dist[14*Np+n];

    		// q=14
    		m14 = dist[13*Np+n];

    		// q=15
    		m15 = dist[16*Np+n];

    		// q=16
    		m16 = dist[15*Np+n];

    		// q=17
    		m17 = dist[18*Np+n];

    		// q=18
    		m18 = dist[17*Np+n];

            //compute fluid velocity
            ux = 3.0/rho0*(m1-m2+m7-m8+m9-m10+m11-m12+m13-m14+0.5*(chem*nx+Fx)/3.0);
            uy = 3.0/rho0*(m3-m4+m7-m8-m9+m10+m15-m16+m17-m18+0.5*(chem*ny+Fy)/3.0);
            uz = 3.0/rho0*(m5-m6+m11-m12-m13+m14+m15-m16-m17+m18+0.5*(chem*nz+Fz)/3.0);

            //compute pressure
            p = (m0+m2+m1+m4+m3+m6+m5+m8+m7+m10+m9+m12+m11+m14+m13+m16+m15+m18+m17)
                      +0.5*(rhoA-rhoB)/2.0/3.0*(ux*nx+uy*ny+uz*nz);

            //compute equilibrium distributions
                feq0 = 0.3333333333333333*p - 0.25*(Fx*ux + Fy*uy + Fz*uz)*(-0.6666666666666666 + ux*ux + uy*uy + uz*uz) - 
         0.16666666666666666*rho0*(ux*ux + uy*uy + uz*uz) - 0.5*(-(nx*ux) - ny*uy - nz*uz)*
          (-0.08333333333333333*(rhoA - rhoB)*(ux*ux + uy*uy + uz*uz) + chem*(0.3333333333333333 - 0.5*(ux*ux + uy*uy + uz*uz)));
                feq1 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-ux*ux + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*uz)) - 
         0.125*(Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 
           0.3333333333333333*(-2.*ux + ux*ux + uy*uy + uz*uz)) - 0.0625*(nx - nx*ux - ny*uy - nz*uz)*
          (2*chem*ux*ux - 0.3333333333333333*((-rhoA + rhoB)*ux*ux + 2*chem*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
           0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*uz)));
                feq2 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-ux*ux + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*uz)) - 
         0.125*(Fx + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 
           0.3333333333333333*(2.*ux + ux*ux + uy*uy + uz*uz)) - 0.0625*(nx + nx*ux + ny*uy + nz*uz)*
          (-2.*chem*ux*ux + 0.1111111111111111*(-4.*chem + rhoB*(-2.*ux - 1.*ux*ux - 1.*uy*uy - 1.*uz*uz) + 
             rhoA*(2.*ux + ux*ux + uy*uy + uz*uz)) + 0.3333333333333333*((-1.*rhoA + rhoB)*ux*ux + 
             chem*(4.*ux + 2.*ux*ux + 2.*uy*uy + 2.*uz*uz)));
                feq3 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uy*uy + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*uz)) - 
         0.125*(Fx*ux + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy + 
           0.3333333333333333*(ux*ux - 2.*uy + uy*uy + uz*uz)) - 0.0625*(ny - nx*ux - ny*uy - nz*uz)*
          (2*chem*uy*uy - 0.3333333333333333*((-rhoA + rhoB)*uy*uy + 2*chem*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
           0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*uz)));
                feq4 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uy*uy + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*uz)) - 
         0.125*(Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy + 
           0.3333333333333333*(ux*ux + 2.*uy + uy*uy + uz*uz)) - 0.0625*(ny + nx*ux + ny*uy + nz*uz)*
          (-2.*chem*uy*uy + 0.1111111111111111*(-4.*chem + rhoB*(-1.*ux*ux - 2.*uy - 1.*uy*uy - 1.*uz*uz) + 
             rhoA*(ux*ux + 2.*uy + uy*uy + uz*uz)) + 0.3333333333333333*((-1.*rhoA + rhoB)*uy*uy + 
             chem*(2.*ux*ux + 4.*uy + 2.*uy*uy + 2.*uz*uz))); 
                feq5 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uz*uz + 0.3333333333333333*(ux*ux + uy*uy + (-2 + uz)*uz)) - 
         0.125*(Fx*ux + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*uz*uz + 
           0.3333333333333333*(ux*ux + uy*uy + (-2. + uz)*uz)) - 0.0625*(nx*ux + ny*uy + nz*(-1. + uz))*
          (-2.*chem*uz*uz + 0.1111111111111111*(-4.*chem + rhoB*(-1.*ux*ux - 1.*uy*uy + (2. - 1.*uz)*uz) + 
             rhoA*(ux*ux + uy*uy + (-2. + uz)*uz)) + 0.3333333333333333*((-1.*rhoA + rhoB)*uz*uz + 
             chem*(2.*ux*ux + 2.*uy*uy + uz*(-4. + 2.*uz)))); 
                feq6 = 0.05555555555555555*p - 0.08333333333333333*rho0*(-uz*uz + 0.3333333333333333*(ux*ux + uy*uy + uz*(2 + uz))) - 
         0.125*(Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*uz*uz + 
           0.3333333333333333*(ux*ux + uy*uy + uz*(2. + uz))) - 0.0625*(nz + nx*ux + ny*uy + nz*uz)*
          (-2.*chem*uz*uz + 0.1111111111111111*(-4.*chem + rhoB*(-1.*ux*ux - 1.*uy*uy + (-2. - 1.*uz)*uz) + 
             rhoA*(ux*ux + uy*uy + uz*(2. + uz))) + 0.3333333333333333*((-1.*rhoA + rhoB)*uz*uz + 
             chem*(2.*ux*ux + 2.*uy*uy + uz*(4. + 2.*uz)))); 
                feq7 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux + uy)*(ux + uy) + 0.3333333333333333*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) - 
         0.0625*(Fx*(-1. + ux) + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uy - 1.*uy*uy + 
           0.3333333333333333*(-2.*ux + ux*ux - 2.*uy + uy*uy + uz*uz)) - 0.03125*(nx + ny - nx*ux - ny*uy - nz*uz)*
          (2*chem*(ux + uy)*(ux + uy) + 0.3333333333333333*((rhoA - rhoB)*(ux + uy)*(ux + uy) - 2*chem*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
           0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz))); 
                feq8 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux + uy)*(ux + uy) + 0.3333333333333333*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) - 
         0.0625*(Fx + Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uy - 1.*uy*uy + 
           0.3333333333333333*(2.*ux + ux*ux + 2.*uy + uy*uy + uz*uz)) - 0.03125*(-(nx*(1 + ux)) - ny*(1 + uy) - nz*uz)*
          (2*chem*(ux + uy)*(ux + uy) - 0.3333333333333333*(-((rhoA - rhoB)*(ux + uy)*(ux + uy)) + 
             2*chem*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz))); 
                feq9 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux - uy)*(ux - uy) + 0.3333333333333333*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) - 
         0.0625*(Fy + Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uy - 1.*uy*uy + 
           0.3333333333333333*(-2.*ux + ux*ux + 2.*uy + uy*uy + uz*uz)) - 0.03125*(nx - nx*ux - ny*(1 + uy) - nz*uz)*
          (2*chem*(ux - uy)*(ux - uy) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uy)*(ux - uy)) + 
             2*chem*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz))); 
                feq10 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux - uy)*(ux - uy) + 0.3333333333333333*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) - 
         0.0625*(Fx*(1 + ux) + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uy - 1.*uy*uy + 
           0.3333333333333333*(2.*ux + ux*ux - 2.*uy + uy*uy + uz*uz)) - 0.03125*(ny - nx*(1 + ux) - ny*uy - nz*uz)*
          (2*chem*(ux - uy)*(ux - uy) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uy)*(ux - uy)) + 
             2*chem*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz))); 
                feq11 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux + uz)*(ux + uz) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) - 
         0.0625*(Fx*(-1. + ux) + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uz - 1.*uz*uz + 
           0.3333333333333333*(-2.*ux + ux*ux + uy*uy + (-2. + uz)*uz)) - 0.03125*(nx + nz - nx*ux - ny*uy - nz*uz)*
          (2*chem*(ux + uz)*(ux + uz) + 0.3333333333333333*((rhoA - rhoB)*(ux + uz)*(ux + uz) - 2*chem*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
           0.1111111111111111*(4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz))); 
                feq12 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux + uz)*(ux + uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) - 
         0.0625*(Fx + Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux - 2.*ux*uz - 1.*uz*uz + 
           0.3333333333333333*(2.*ux + ux*ux + uy*uy + uz*(2. + uz))) - 0.03125*(-(nx*(1 + ux)) - ny*uy - nz*(1 + uz))*
          (2*chem*(ux + uz)*(ux + uz) - 0.3333333333333333*(-((rhoA - rhoB)*(ux + uz)*(ux + uz)) + 
             2*chem*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))); 
                feq13 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux - uz)*(ux - uz) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) - 
         0.0625*(Fz + Fx*(-1. + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uz - 1.*uz*uz + 
           0.3333333333333333*(-2.*ux + ux*ux + uy*uy + uz*(2. + uz))) - 0.03125*(nx - nx*ux - ny*uy - nz*(1 + uz))*
          (2*chem*(ux - uz)*(ux - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uz)*(ux - uz)) + 
             2*chem*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*(2 + uz)))); 
                feq14 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(ux - uz)*(ux - uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) - 
         0.0625*(Fx*(1 + ux) + Fy*uy + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*ux*ux + 2.*ux*uz - 1.*uz*uz + 
           0.3333333333333333*(2.*ux + ux*ux + uy*uy + (-2. + uz)*uz)) - 0.03125*(nz - nx*(1 + ux) - ny*uy - nz*uz)*
          (2*chem*(ux - uz)*(ux - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(ux - uz)*(ux - uz)) + 
             2*chem*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz))); 
                feq15 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) - 
         0.0625*(Fx*ux + Fy*(-1. + uy) + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*uy*uy - 2.*uy*uz - 1.*uz*uz + 
           0.3333333333333333*(ux*ux - 2.*uy + uy*uy + (-2. + uz)*uz)) - 0.03125*(ny + nz - nx*ux - ny*uy - nz*uz)*
          (2*chem*(uy + uz)*(uy + uz) + 0.3333333333333333*((rhoA - rhoB)*(uy + uz)*(uy + uz) - 2*chem*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) + 
           0.1111111111111111*(4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))); 
                feq16 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) - 
         0.0625*(Fy + Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy - 2.*uy*uz - 1.*uz*uz + 
           0.3333333333333333*(ux*ux + 2.*uy + uy*uy + uz*(2. + uz))) - 0.03125*(-(nx*ux) - ny*(1 + uy) - nz*(1 + uz))*
          (2*chem*(uy + uz)*(uy + uz) - 0.3333333333333333*(-((rhoA - rhoB)*(uy + uz)*(uy + uz)) + 
             2*chem*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))); 
                feq17 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) - 
         0.0625*(Fz + Fx*ux + Fy*(-1. + uy) + Fz*uz)*(-0.2222222222222222 - 1.*uy*uy + 2.*uy*uz - 1.*uz*uz + 
           0.3333333333333333*(ux*ux - 2.*uy + uy*uy + uz*(2. + uz))) - 0.03125*(ny - nx*ux - ny*uy - nz*(1 + uz))*
          (2*chem*(uy - uz)*(uy - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(uy - uz)*(uy - uz)) + 
             2*chem*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))); 
                feq18 = 0.027777777777777776*p - 0.041666666666666664*rho0*
          (-(uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) - 
         0.0625*(Fx*ux + Fy*(1 + uy) + Fz*(-1. + uz))*(-0.2222222222222222 - 1.*uy*uy + 2.*uy*uz - 1.*uz*uz + 
           0.3333333333333333*(ux*ux + 2.*uy + uy*uy + (-2. + uz)*uz)) - 0.03125*(nz - nx*ux - ny*(1 + uy) - nz*uz)*
          (2*chem*(uy - uz)*(uy - uz) - 0.3333333333333333*(-((rhoA - rhoB)*(uy - uz)*(uy - uz)) + 
             2*chem*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 0.1111111111111111*
            (4*chem - (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz))); 

            //------------------------------------------------- BCK collison ------------------------------------------------------------//
    		// q=0
    		dist[n] = m0 - (m0-feq0)/tau + 0.25*(2*(Fx*ux + Fy*uy + Fz*uz)*(-0.6666666666666666 + ux*ux + uy*uy + uz*uz) + 
         (nx*ux + ny*uy + nz*uz)*(2*chem*(ux*ux + uy*uy + uz*uz) + 
            0.3333333333333333*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + uz*uz))));

    		// q = 1
    		dist[1*Np+n] = m1 - (m1-feq1)/tau + 0.125*(2*(Fx*(-1 + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - ux*ux + 
            0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
         (nx*(-1 + ux) + ny*uy + nz*uz)*(-2*chem*(ux*ux) + 
            0.3333333333333333*((-rhoA + rhoB)*(ux*ux) + 2*chem*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*uz))));

    		// q=2
    		dist[2*Np+n] = m2 - (m2-feq2)/tau + 0.125*(2*(Fx + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - ux*ux + 
            0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*uz)) + 
         (nx + nx*ux + ny*uy + nz*uz)*(-2*chem*(ux*ux) + 
            0.3333333333333333*((-rhoA + rhoB)*(ux*ux) + 2*chem*(2*ux + ux*ux + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*uz))));

    		// q = 3
    		dist[3*Np+n] = m3 - (m3-feq3)/tau + 0.125*(2*(Fx*ux + Fy*(-1 + uy) + Fz*uz)*(-0.2222222222222222 - uy*uy + 
            0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
         (nx*ux + ny*(-1 + uy) + nz*uz)*(-2*chem*(uy*uy) + 
            0.3333333333333333*((-rhoA + rhoB)*(uy*uy) + 2*chem*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*uz))));

    		// q = 4
    		dist[4*Np+n] = m4 - (m4-feq4)/tau + 0.125*(2*(Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - uy*uy + 
            0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*uz)) + 
         (ny + nx*ux + ny*uy + nz*uz)*(-2*chem*(uy*uy) + 
            0.3333333333333333*((-rhoA + rhoB)*(uy*uy) + 2*chem*(ux*ux + 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*uz))));

    		// q = 5
    		dist[5*Np+n] = m5 - (m5-feq5)/tau + 0.125*(2*(Fx*ux + Fy*uy + Fz*(-1 + uz))*(-0.2222222222222222 - uz*uz + 
            0.3333333333333333*(ux*ux + uy*uy + (-2 + uz)*uz)) + 
         (nx*ux + ny*uy + nz*(-1 + uz))*(-2*chem*(uz*uz) + 
            0.3333333333333333*((-rhoA + rhoB)*(uz*uz) + 2*chem*(ux*ux + uy*uy + (-2 + uz)*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + (-2 + uz)*uz))));

    		// q = 6
    		dist[6*Np+n] = m6 - (m6-feq6)/tau + 0.125*(2*(Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - uz*uz + 
            0.3333333333333333*(ux*ux + uy*uy + uz*(2 + uz))) + 
         (nz + nx*ux + ny*uy + nz*uz)*(-2*chem*(uz*uz) + 
            0.3333333333333333*((-rhoA + rhoB)*(uz*uz) + 2*chem*(ux*ux + uy*uy + uz*(2 + uz))) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + uy*uy + uz*(2 + uz)))));

    		// q = 7
    		dist[7*Np+n] = m7 - (m7-feq7)/tau + 0.0625*(-2*(Fx*(-1 + ux) + Fy*(-1 + uy) + Fz*uz)*
          (0.2222222222222222 + (ux + uy)*(ux + uy) - 
            0.3333333333333333*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
         (nx*(-1 + ux) + ny*(-1 + uy) + nz*uz)*
          (-2*chem*((ux + uy)*(ux + uy)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux + uy)*(ux + uy))) + 2*chem*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz))));

    		// q = 8
    		dist[8*Np+n] = m8 - (m8-feq8)/tau + 0.0625*(2*(Fx + Fy + Fx*ux + Fy*uy + Fz*uz)*
          (-0.2222222222222222 - (ux + uy)*(ux + uy) + 
            0.3333333333333333*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
         (nx + ny + nx*ux + ny*uy + nz*uz)*
          (-2*chem*((ux + uy)*(ux + uy)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux + uy)*(ux + uy))) + 2*chem*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz))));

    		// q = 9
    		dist[9*Np+n] = m9 - (m9-feq9)/tau + 0.0625*(2*(Fy + Fx*(-1 + ux) + Fy*uy + Fz*uz)*
          (-0.2222222222222222 - (ux - uy)*(ux - uy) + 
            0.3333333333333333*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
         (ny + nx*(-1 + ux) + ny*uy + nz*uz)*
          (-2*chem*((ux - uy)*(ux - uy)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux - uy)*(ux - uy))) + 2*chem*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz))));

    		// q = 10
    		dist[10*Np+n] = m10 - (m10-feq10)/tau + 0.0625*(2*(Fx*(1 + ux) + Fy*(-1 + uy) + Fz*uz)*
          (-0.2222222222222222 - (ux - uy)*(ux - uy) + 
            0.3333333333333333*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
         (nx*(1 + ux) + ny*(-1 + uy) + nz*uz)*
          (-2*chem*((ux - uy)*(ux - uy)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux - uy)*(ux - uy))) + 2*chem*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz))));

    		// q = 11
    		dist[11*Np+n] = m11 - (m11-feq11)/tau + 0.0625*(-2*(Fx*(-1 + ux) + Fy*uy + Fz*(-1 + uz))*
          (0.2222222222222222 + (ux + uz)*(ux + uz) - 
            0.3333333333333333*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
         (nx*(-1 + ux) + ny*uy + nz*(-1 + uz))*
          (-2*chem*((ux + uz)*(ux + uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux + uz)*(ux + uz))) + 2*chem*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz))));

    		// q = 12
    		dist[12*Np+n] = m12 - (m12-feq12)/tau + 0.0625*(2*(Fx + Fz + Fx*ux + Fy*uy + Fz*uz)*
          (-0.2222222222222222 - (ux + uz)*(ux + uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))
           + (nx + nz + nx*ux + ny*uy + nz*uz)*
          (-2*chem*((ux + uz)*(ux + uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux + uz)*(ux + uz))) + 2*chem*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))));

    		// q = 13
    		dist[13*Np+n] = m13 - (m13-feq13)/tau + 0.0625*(2*(Fz + Fx*(-1 + ux) + Fy*uy + Fz*uz)*
          (-0.2222222222222222 - (ux - uz)*(ux - uz) + 
            0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
         (nz + nx*(-1 + ux) + ny*uy + nz*uz)*
          (-2*chem*((ux - uz)*(ux - uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux - uz)*(ux - uz))) + 2*chem*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(-2*ux + ux*ux + uy*uy + uz*(2 + uz)))));

    		// q= 14
    		dist[14*Np+n] = m14 - (m14-feq14)/tau + 0.0625*(2*(Fx*(1 + ux) + Fy*uy + Fz*(-1 + uz))*
          (-0.2222222222222222 - (ux - uz)*(ux - uz) + 
            0.3333333333333333*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
         (nx*(1 + ux) + ny*uy + nz*(-1 + uz))*
          (-2*chem*((ux - uz)*(ux - uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((ux - uz)*(ux - uz))) + 2*chem*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz))));

    		// q = 15
    		dist[15*Np+n] = m15 - (m15-feq15)/tau + 0.0625*(-2*(Fx*ux + Fy*(-1 + uy) + Fz*(-1 + uz))*
          (0.2222222222222222 + (uy + uz)*(uy + uz) - 0.3333333333333333*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))
           + (nx*ux + ny*(-1 + uy) + nz*(-1 + uz))*
          (-2*chem*((uy + uz)*(uy + uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((uy + uz)*(uy + uz))) + 2*chem*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))));

    		// q = 16
    		dist[16*Np+n] = m16 - (m16-feq16)/tau + 0.0625*(2*(Fy + Fz + Fx*ux + Fy*uy + Fz*uz)*
          (-0.2222222222222222 - (uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))
           + (ny + nz + nx*ux + ny*uy + nz*uz)*
          (-2*chem*((uy + uz)*(uy + uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((uy + uz)*(uy + uz))) + 2*chem*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))));

    		// q = 17
    		dist[17*Np+n] = m17 - (m17-feq17)/tau + 0.0625*(2*(Fz + Fx*ux + Fy*(-1 + uy) + Fz*uz)*
          (-0.2222222222222222 - (uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))
           + (nz + nx*ux + ny*(-1 + uy) + nz*uz)*
          (-2*chem*((uy - uz)*(uy - uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((uy - uz)*(uy - uz))) + 2*chem*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))));

    		// q = 18
    		dist[18*Np+n] = m18 - (m18-feq18)/tau + 0.0625*(2*(Fx*ux + Fy*(1 + uy) + Fz*(-1 + uz))*
          (-0.2222222222222222 - (uy - uz)*(uy - uz) + 
            0.3333333333333333*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 
         (nx*ux + ny*(1 + uy) + nz*(-1 + uz))*
          (-2*chem*((uy - uz)*(uy - uz)) + 0.3333333333333333*
             (-((rhoA - rhoB)*((uy - uz)*(uy - uz))) + 2*chem*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 
            0.1111111111111111*(-4*chem + (rhoA - rhoB)*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz))));
            //----------------------------------------------------------------------------------------------------------------------------------------//

            // ----------------------------- compute phase field evolution ----------------------------------------
            //Normalize the Color Gradient
            C = sqrt(nx*nx+ny*ny+nz*nz);
            double ColorMag = C;
            if (C==0.0) ColorMag=1.0;
            nx = nx/ColorMag;
            ny = ny/ColorMag;
            nz = nz/ColorMag;		
            //compute surface tension-related parameter
            //theta = 4.5*M*2.0*(1-phi*phi)/W;
            theta = 4.5*M*2.0*(1-phi_temp*phi_temp)/W;

            //load distributions of phase field
            //q=0
            h0 = hq[n];
            //q=1
            h1 = hq[2*Np+n]; 

            //q=2
            h2 = hq[1*Np+n];  

            //q=3
            h3 = hq[4*Np+n];

            //q=4
            h4 = hq[3*Np+n];

            //q=5
            h5 = hq[6*Np+n];

            //q=6
            h6 = hq[5*Np+n];

            //-------------------------------- BGK collison for phase field ---------------------------------//
            // q = 0
            hq[n] = h0 - (h0 - 0.3333333333333333*phi)/tauM;

            // q = 1
            hq[1*Np+n] = h1 - (h1 - 0.1111111111111111*nx*theta - phi*(0.1111111111111111 + 0.5*ux))/tauM;

            // q = 2
            hq[2*Np+n] = h2 - (h2 + 0.1111111111111111*nx*theta - phi*(0.1111111111111111 - 0.5*ux))/tauM;

            // q = 3
            hq[3*Np+n] = h3 - (h3 - 0.1111111111111111*ny*theta - phi*(0.1111111111111111 + 0.5*uy))/tauM;

            // q = 4
            hq[4*Np+n] = h4 - (h4 + 0.1111111111111111*ny*theta - phi*(0.1111111111111111 - 0.5*uy))/tauM;

            // q = 5
            hq[5*Np+n] = h5 - (h5 - 0.1111111111111111*nz*theta - phi*(0.1111111111111111 + 0.5*uz))/tauM;

            // q = 6
            hq[6*Np+n] = h6 - (h6 + 0.1111111111111111*nz*theta - phi*(0.1111111111111111 - 0.5*uz))/tauM;
            //........................................................................

            //Update velocity on device
    		Vel[0*Np+n] = ux;
    		Vel[1*Np+n] = uy;
    		Vel[2*Np+n] = uz;
            //Update pressure on device
            Pressure[n] = p;
            //Update chemical potential on device
            mu_phi[n] = chem;
            //Update color gradient on device
    		ColorGrad[0*Np+n] = nx;
    		ColorGrad[1*Np+n] = ny;
    		ColorGrad[2*Np+n] = nz;

    	}
	}
}

__global__ void dvc_ScaLBL_D3Q19_AAodd_FreeLeeModel_SingleFluid_BGK(int *neighborList, double *dist, double *Vel, double *Pressure,  
		double tau, double rho0, double Fx, double Fy, double Fz, int start, int finish, int Np){

	int n;
	int nr1,nr2,nr3,nr4,nr5,nr6,nr7,nr8,nr9,nr10,nr11,nr12,nr13,nr14,nr15,nr16,nr17,nr18;
	double ux,uy,uz;//fluid velocity 
	double p;//pressure
	// distribution functions
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
	double m0,m3,m5,m7;

	//	for (int n=start; n<finish; n++){
	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

		if ( n<finish ){
			// q=0
			m0 = dist[n];
			// q=1
			nr1 = neighborList[n]; // neighbor 2 ( > 10Np => odd part of dist)
			m1 = dist[nr1]; // reading the f1 data into register fq

			nr2 = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
			m2 = dist[nr2];  // reading the f2 data into register fq

			// q=3
			nr3 = neighborList[n+2*Np]; // neighbor 4
			m3 = dist[nr3];

			// q = 4
			nr4 = neighborList[n+3*Np]; // neighbor 3
			m4 = dist[nr4];

			// q=5
			nr5 = neighborList[n+4*Np];
			m5 = dist[nr5];

			// q = 6
			nr6 = neighborList[n+5*Np];
			m6 = dist[nr6];

			// q=7
			nr7 = neighborList[n+6*Np];
			m7 = dist[nr7];

			// q = 8
			nr8 = neighborList[n+7*Np];
			m8 = dist[nr8];

			// q=9
			nr9 = neighborList[n+8*Np];
			m9 = dist[nr9];

			// q = 10
			nr10 = neighborList[n+9*Np];
			m10 = dist[nr10];

			// q=11
			nr11 = neighborList[n+10*Np];
			m11 = dist[nr11];

			// q=12
			nr12 = neighborList[n+11*Np];
			m12 = dist[nr12];

			// q=13
			nr13 = neighborList[n+12*Np];
			m13 = dist[nr13];

			// q=14
			nr14 = neighborList[n+13*Np];
			m14 = dist[nr14];

			// q=15
			nr15 = neighborList[n+14*Np];
			m15 = dist[nr15];

			// q=16
			nr16 = neighborList[n+15*Np];
			m16 = dist[nr16];

			// q=17
			nr17 = neighborList[n+16*Np];
			m17 = dist[nr17];

			// q=18
			nr18 = neighborList[n+17*Np];
			m18 = dist[nr18];

			//compute fluid velocity
			ux = 3.0/rho0*(m1-m2+m7-m8+m9-m10+m11-m12+m13-m14+0.5*(Fx)/3.0);
			uy = 3.0/rho0*(m3-m4+m7-m8-m9+m10+m15-m16+m17-m18+0.5*(Fy)/3.0);
			uz = 3.0/rho0*(m5-m6+m11-m12-m13+m14+m15-m16-m17+m18+0.5*(Fz)/3.0);
			//compute pressure
			p = (m0+m2+m1+m4+m3+m6+m5+m8+m7+m10+m9+m12+m11+m14+m13+m16+m15+m18+m17);

			//------------------------------------------------- BCK collison ------------------------------------------------------------//
			// q=0
			dist[n] = m0 + 0.5*(Fx*ux + Fy*uy + Fz*uz)*(-0.6666666666666666 + ux*ux + uy*uy + uz*uz) - 
					(m0 - 0.3333333333333333*p + 0.25*(Fx*ux + Fy*uy + Fz*uz)*
							(-0.6666666666666666 + ux*ux + uy*uy + uz*uz) + 0.16666666666666666*rho0*(ux*ux + uy*uy + uz*uz))/
							tau;

			// q = 1
			dist[nr2] = m1 + 0.25*(Fx*(-1 + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - ux*ux + 
					0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*uz)) - 
					(m1 - 0.05555555555555555*p + 0.08333333333333333*rho0*
							(-(ux*ux) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
							0.125*(Fx*(-1. + ux) + Fy*uy + Fz*uz)*
							(-0.2222222222222222 - 1.*(ux*ux) + 0.3333333333333333*(-2.*ux + ux*ux + uy*uy + uz*uz)))/tau;

			// q=2
			dist[nr1] = m2 + 0.25*(Fx + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - ux*ux + 
					0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*uz)) - 
					(m2 - 0.05555555555555555*p + 0.08333333333333333*rho0*
							(-(ux*ux) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*uz)) + 
							0.125*(Fx + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(ux*ux) + 
									0.3333333333333333*(2.*ux + ux*ux + uy*uy + uz*uz)))/tau;

			// q = 3
			dist[nr4] = m3 + 0.25*(Fx*ux + Fy*(-1 + uy) + Fz*uz)*(-0.2222222222222222 - uy*uy + 
					0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*uz)) - 
					(m3 - 0.05555555555555555*p + 0.08333333333333333*rho0*
							(-(uy*uy) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
							0.125*(Fx*ux + Fy*(-1. + uy) + Fz*uz)*
							(-0.2222222222222222 - 1.*(uy*uy) + 0.3333333333333333*(ux*ux - 2.*uy + uy*uy + uz*uz)))/tau;

			// q = 4
			dist[nr3] = m4 + 0.25*(Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - uy*uy + 
					0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*uz)) - 
					(m4 - 0.05555555555555555*p + 0.08333333333333333*rho0*
							(-(uy*uy) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*uz)) + 
							0.125*(Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(uy*uy) + 
									0.3333333333333333*(ux*ux + 2.*uy + uy*uy + uz*uz)))/tau;

			// q = 5
			dist[nr6] = m5 + 0.25*(Fx*ux + Fy*uy + Fz*(-1 + uz))*(-0.2222222222222222 - uz*uz + 
					0.3333333333333333*(ux*ux + uy*uy + (-2 + uz)*uz)) - 
					(m5 - 0.05555555555555555*p + 0.08333333333333333*rho0*
							(-(uz*uz) + 0.3333333333333333*(ux*ux + uy*uy + (-2 + uz)*uz)) + 
							0.125*(Fx*ux + Fy*uy + Fz*(-1. + uz))*
							(-0.2222222222222222 - 1.*(uz*uz) + 0.3333333333333333*(ux*ux + uy*uy + (-2. + uz)*uz)))/tau;

			// q = 6
			dist[nr5] = m6 + 0.25*(Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - uz*uz + 
					0.3333333333333333*(ux*ux + uy*uy + uz*(2 + uz))) - 
					(m6 - 0.05555555555555555*p + 0.08333333333333333*rho0*
							(-(uz*uz) + 0.3333333333333333*(ux*ux + uy*uy + uz*(2 + uz))) + 
							0.125*(Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(uz*uz) + 
									0.3333333333333333*(ux*ux + uy*uy + uz*(2. + uz))))/tau;

			// q = 7
			dist[nr8] = m7 - 0.125*(Fx*(-1 + ux) + Fy*(-1 + uy) + Fz*uz)*
					(0.2222222222222222 + (ux + uy)*(ux + uy) - 0.3333333333333333*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz))\
					- (m7 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((ux + uy)*(ux + uy)) + 0.3333333333333333*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
							0.0625*(Fx*(-1. + ux) + Fy*(-1. + uy) + Fz*uz)*
							(-0.2222222222222222 - 1.*(ux*ux) - 2.*ux*uy - 1.*(uy*uy) + 
									0.3333333333333333*(-2.*ux + ux*ux - 2.*uy + uy*uy + uz*uz)))/tau;

			// q = 8
			dist[nr7] = m8 + 0.125*(Fx + Fy + Fx*ux + Fy*uy + Fz*uz)*
					(-0.2222222222222222 - (ux + uy)*(ux + uy) + 0.3333333333333333*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz))\
					- (m8 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((ux + uy)*(ux + uy)) + 0.3333333333333333*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
							0.0625*(Fx + Fy + Fx*ux + Fy*uy + Fz*uz)*
							(-0.2222222222222222 - 1.*(ux*ux) - 2.*ux*uy - 1.*(uy*uy) + 
									0.3333333333333333*(2.*ux + ux*ux + 2.*uy + uy*uy + uz*uz)))/tau;

			// q = 9
			dist[nr10] = m9 + 0.125*(Fy + Fx*(-1 + ux) + Fy*uy + Fz*uz)*
					(-0.2222222222222222 - (ux - uy)*(ux - uy) + 0.3333333333333333*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz))
					- (m9 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((ux - uy)*(ux - uy)) + 0.3333333333333333*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
							0.0625*(Fy + Fx*(-1. + ux) + Fy*uy + Fz*uz)*
							(-0.2222222222222222 - 1.*(ux*ux) + 2.*ux*uy - 1.*(uy*uy) + 
									0.3333333333333333*(-2.*ux + ux*ux + 2.*uy + uy*uy + uz*uz)))/tau;

			// q = 10
			dist[nr9] = m10 + 0.125*(Fx*(1 + ux) + Fy*(-1 + uy) + Fz*uz)*
					(-0.2222222222222222 - (ux - uy)*(ux - uy) + 0.3333333333333333*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz))\
					- (m10 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((ux - uy)*(ux - uy)) + 0.3333333333333333*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
							0.0625*(Fx*(1 + ux) + Fy*(-1. + uy) + Fz*uz)*
							(-0.2222222222222222 - 1.*(ux*ux) + 2.*ux*uy - 1.*(uy*uy) + 
									0.3333333333333333*(2.*ux + ux*ux - 2.*uy + uy*uy + uz*uz)))/tau;

			// q = 11
			dist[nr12] = m11 - 0.125*(Fx*(-1 + ux) + Fy*uy + Fz*(-1 + uz))*
					(0.2222222222222222 + (ux + uz)*(ux + uz) - 0.3333333333333333*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz))\
					- (m11 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((ux + uz)*(ux + uz)) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
							0.0625*(Fx*(-1. + ux) + Fy*uy + Fz*(-1. + uz))*
							(-0.2222222222222222 - 1.*(ux*ux) - 2.*ux*uz - 1.*(uz*uz) + 
									0.3333333333333333*(-2.*ux + ux*ux + uy*uy + (-2. + uz)*uz)))/tau;

			// q = 12
			dist[nr11] = m12 + 0.125*(Fx + Fz + Fx*ux + Fy*uy + Fz*uz)*
					(-0.2222222222222222 - (ux + uz)*(ux + uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))\
					- (m12 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((ux + uz)*(ux + uz)) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
							0.0625*(Fx + Fz + Fx*ux + Fy*uy + Fz*uz)*
							(-0.2222222222222222 - 1.*(ux*ux) - 2.*ux*uz - 1.*(uz*uz) + 
									0.3333333333333333*(2.*ux + ux*ux + uy*uy + uz*(2. + uz))))/tau;

			// q = 13
			dist[nr14] = m13 + 0.125*(Fz + Fx*(-1 + ux) + Fy*uy + Fz*uz)*
					(-0.2222222222222222 - (ux - uz)*(ux - uz) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*(2 + uz)))\
					- (m13 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((ux - uz)*(ux - uz)) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
							0.0625*(Fz + Fx*(-1. + ux) + Fy*uy + Fz*uz)*
							(-0.2222222222222222 - 1.*(ux*ux) + 2.*ux*uz - 1.*(uz*uz) + 
									0.3333333333333333*(-2.*ux + ux*ux + uy*uy + uz*(2. + uz))))/tau;

			// q= 14
			dist[nr13] = m14 + 0.125*(Fx*(1 + ux) + Fy*uy + Fz*(-1 + uz))*
					(-0.2222222222222222 - (ux - uz)*(ux - uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz))\
					- (m14 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((ux - uz)*(ux - uz)) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
							0.0625*(Fx*(1 + ux) + Fy*uy + Fz*(-1. + uz))*
							(-0.2222222222222222 - 1.*(ux*ux) + 2.*ux*uz - 1.*(uz*uz) + 
									0.3333333333333333*(2.*ux + ux*ux + uy*uy + (-2. + uz)*uz)))/tau;

			// q = 15
			dist[nr16] = m15 - 0.125*(Fx*ux + Fy*(-1 + uy) + Fz*(-1 + uz))*
					(0.2222222222222222 + (uy + uz)*(uy + uz) - 0.3333333333333333*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))\
					- (m15 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((uy + uz)*(uy + uz)) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) + 
							0.0625*(Fx*ux + Fy*(-1. + uy) + Fz*(-1. + uz))*
							(-0.2222222222222222 - 1.*(uy*uy) - 2.*uy*uz - 1.*(uz*uz) + 
									0.3333333333333333*(ux*ux - 2.*uy + uy*uy + (-2. + uz)*uz)))/tau;

			// q = 16
			dist[nr15] = m16 + 0.125*(Fy + Fz + Fx*ux + Fy*uy + Fz*uz)*
					(-0.2222222222222222 - (uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))\
					- (m16 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((uy + uz)*(uy + uz)) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) + 
							0.0625*(Fy + Fz + Fx*ux + Fy*uy + Fz*uz)*
							(-0.2222222222222222 - 1.*(uy*uy) - 2.*uy*uz - 1.*(uz*uz) + 
									0.3333333333333333*(ux*ux + 2.*uy + uy*uy + uz*(2. + uz))))/tau;

			// q = 17
			dist[nr18] = m17 + 0.125*(Fz + Fx*ux + Fy*(-1 + uy) + Fz*uz)*
					(-0.2222222222222222 - (uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))\
					- (m17 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((uy - uz)*(uy - uz)) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) + 
							0.0625*(Fz + Fx*ux + Fy*(-1. + uy) + Fz*uz)*
							(-0.2222222222222222 - 1.*(uy*uy) + 2.*uy*uz - 1.*(uz*uz) + 
									0.3333333333333333*(ux*ux - 2.*uy + uy*uy + uz*(2. + uz))))/tau;

			// q = 18
			dist[nr17] = m18 + 0.125*(Fx*ux + Fy*(1 + uy) + Fz*(-1 + uz))*
					(-0.2222222222222222 - (uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz))\
					- (m18 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((uy - uz)*(uy - uz)) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 
							0.0625*(Fx*ux + Fy*(1 + uy) + Fz*(-1. + uz))*
							(-0.2222222222222222 - 1.*(uy*uy) + 2.*uy*uz - 1.*(uz*uz) + 
									0.3333333333333333*(ux*ux + 2.*uy + uy*uy + (-2. + uz)*uz)))/tau;
			//----------------------------------------------------------------------------------------------------------------------------------------//


			//Update velocity on device
			Vel[0*Np+n] = ux;
			Vel[1*Np+n] = uy;
			Vel[2*Np+n] = uz;
			//Update pressure on device
			Pressure[n] = p;
		}
	}
}

__global__ void dvc_ScaLBL_D3Q19_AAeven_FreeLeeModel_SingleFluid_BGK(double *dist, double *Vel, double *Pressure, 
		double tau, double rho0, double Fx, double Fy, double Fz, int start, int finish, int Np){

	int n;
	double ux,uy,uz;//fluid velocity 
	double p;//pressure
	// distribution functions
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
	double m0,m3,m5,m7;

	//	for (int n=start; n<finish; n++){
	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

		if ( n<finish ){	
			// q=0
			m0 = dist[n];
			// q=1
			m1 = dist[2*Np+n]; 

			// q=2
			m2 = dist[1*Np+n];  

			// q=3
			m3 = dist[4*Np+n];

			// q = 4
			m4 = dist[3*Np+n];

			// q=5
			m5 = dist[6*Np+n];

			// q = 6
			m6 = dist[5*Np+n];

			// q=7
			m7 = dist[8*Np+n];

			// q = 8
			m8 = dist[7*Np+n];

			// q=9
			m9 = dist[10*Np+n];

			// q = 10
			m10 = dist[9*Np+n];

			// q=11
			m11 = dist[12*Np+n];

			// q=12
			m12 = dist[11*Np+n];

			// q=13
			m13 = dist[14*Np+n];

			// q=14
			m14 = dist[13*Np+n];

			// q=15
			m15 = dist[16*Np+n];

			// q=16
			m16 = dist[15*Np+n];

			// q=17
			m17 = dist[18*Np+n];

			// q=18
			m18 = dist[17*Np+n];

			//compute fluid velocity
			ux = 3.0/rho0*(m1-m2+m7-m8+m9-m10+m11-m12+m13-m14+0.5*(Fx)/3.0);
			uy = 3.0/rho0*(m3-m4+m7-m8-m9+m10+m15-m16+m17-m18+0.5*(Fy)/3.0);
			uz = 3.0/rho0*(m5-m6+m11-m12-m13+m14+m15-m16-m17+m18+0.5*(Fz)/3.0);
			//compute pressure
			p = (m0+m2+m1+m4+m3+m6+m5+m8+m7+m10+m9+m12+m11+m14+m13+m16+m15+m18+m17);

			//------------------------------------------------- BCK collison ------------------------------------------------------------//
			// q=0
			dist[n] = m0 + 0.5*(Fx*ux + Fy*uy + Fz*uz)*(-0.6666666666666666 + ux*ux + uy*uy + uz*uz) - 
					(m0 - 0.3333333333333333*p + 0.25*(Fx*ux + Fy*uy + Fz*uz)*
							(-0.6666666666666666 + ux*ux + uy*uy + uz*uz) + 0.16666666666666666*rho0*(ux*ux + uy*uy + uz*uz))/
							tau; 

			// q = 1
			dist[1*Np+n] = m1 + 0.25*(Fx*(-1 + ux) + Fy*uy + Fz*uz)*(-0.2222222222222222 - ux*ux + 
					0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*uz)) - 
					(m1 - 0.05555555555555555*p + 0.08333333333333333*rho0*
							(-(ux*ux) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*uz)) + 
							0.125*(Fx*(-1. + ux) + Fy*uy + Fz*uz)*
							(-0.2222222222222222 - 1.*(ux*ux) + 0.3333333333333333*(-2.*ux + ux*ux + uy*uy + uz*uz)))/tau;

			// q=2
			dist[2*Np+n] = m2 + 0.25*(Fx + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - ux*ux + 
					0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*uz)) - 
					(m2 - 0.05555555555555555*p + 0.08333333333333333*rho0*
							(-(ux*ux) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*uz)) + 
							0.125*(Fx + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(ux*ux) + 
									0.3333333333333333*(2.*ux + ux*ux + uy*uy + uz*uz)))/tau;

			// q = 3
			dist[3*Np+n] = m3 + 0.25*(Fx*ux + Fy*(-1 + uy) + Fz*uz)*(-0.2222222222222222 - uy*uy + 
					0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*uz)) - 
					(m3 - 0.05555555555555555*p + 0.08333333333333333*rho0*
							(-(uy*uy) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*uz)) + 
							0.125*(Fx*ux + Fy*(-1. + uy) + Fz*uz)*
							(-0.2222222222222222 - 1.*(uy*uy) + 0.3333333333333333*(ux*ux - 2.*uy + uy*uy + uz*uz)))/tau;

			// q = 4
			dist[4*Np+n] = m4 + 0.25*(Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - uy*uy + 
					0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*uz)) - 
					(m4 - 0.05555555555555555*p + 0.08333333333333333*rho0*
							(-(uy*uy) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*uz)) + 
							0.125*(Fy + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(uy*uy) + 
									0.3333333333333333*(ux*ux + 2.*uy + uy*uy + uz*uz)))/tau;

			// q = 5
			dist[5*Np+n] = m5 + 0.25*(Fx*ux + Fy*uy + Fz*(-1 + uz))*(-0.2222222222222222 - uz*uz + 
					0.3333333333333333*(ux*ux + uy*uy + (-2 + uz)*uz)) - 
					(m5 - 0.05555555555555555*p + 0.08333333333333333*rho0*
							(-(uz*uz) + 0.3333333333333333*(ux*ux + uy*uy + (-2 + uz)*uz)) + 
							0.125*(Fx*ux + Fy*uy + Fz*(-1. + uz))*
							(-0.2222222222222222 - 1.*(uz*uz) + 0.3333333333333333*(ux*ux + uy*uy + (-2. + uz)*uz)))/tau;

			// q = 6
			dist[6*Np+n] = m6 + 0.25*(Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - uz*uz + 
					0.3333333333333333*(ux*ux + uy*uy + uz*(2 + uz))) - 
					(m6 - 0.05555555555555555*p + 0.08333333333333333*rho0*
							(-(uz*uz) + 0.3333333333333333*(ux*ux + uy*uy + uz*(2 + uz))) + 
							0.125*(Fz + Fx*ux + Fy*uy + Fz*uz)*(-0.2222222222222222 - 1.*(uz*uz) + 
									0.3333333333333333*(ux*ux + uy*uy + uz*(2. + uz))))/tau;

			// q = 7
			dist[7*Np+n] = m7 - 0.125*(Fx*(-1 + ux) + Fy*(-1 + uy) + Fz*uz)*
					(0.2222222222222222 + (ux + uy)*(ux + uy) - 0.3333333333333333*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz))\
					- (m7 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((ux + uy)*(ux + uy)) + 0.3333333333333333*(-2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
							0.0625*(Fx*(-1. + ux) + Fy*(-1. + uy) + Fz*uz)*
							(-0.2222222222222222 - 1.*(ux*ux) - 2.*ux*uy - 1.*(uy*uy) + 
									0.3333333333333333*(-2.*ux + ux*ux - 2.*uy + uy*uy + uz*uz)))/tau;

			// q = 8
			dist[8*Np+n] = m8 + 0.125*(Fx + Fy + Fx*ux + Fy*uy + Fz*uz)*
					(-0.2222222222222222 - (ux + uy)*(ux + uy) + 0.3333333333333333*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz))\
					- (m8 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((ux + uy)*(ux + uy)) + 0.3333333333333333*(2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
							0.0625*(Fx + Fy + Fx*ux + Fy*uy + Fz*uz)*
							(-0.2222222222222222 - 1.*(ux*ux) - 2.*ux*uy - 1.*(uy*uy) + 
									0.3333333333333333*(2.*ux + ux*ux + 2.*uy + uy*uy + uz*uz)))/tau;

			// q = 9
			dist[9*Np+n] = m9 + 0.125*(Fy + Fx*(-1 + ux) + Fy*uy + Fz*uz)*
					(-0.2222222222222222 - (ux - uy)*(ux - uy) + 0.3333333333333333*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz))
					- (m9 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((ux - uy)*(ux - uy)) + 0.3333333333333333*(-2*ux + ux*ux + 2*uy + uy*uy + uz*uz)) + 
							0.0625*(Fy + Fx*(-1. + ux) + Fy*uy + Fz*uz)*
							(-0.2222222222222222 - 1.*(ux*ux) + 2.*ux*uy - 1.*(uy*uy) + 
									0.3333333333333333*(-2.*ux + ux*ux + 2.*uy + uy*uy + uz*uz)))/tau;

			// q = 10
			dist[10*Np+n] = m10 + 0.125*(Fx*(1 + ux) + Fy*(-1 + uy) + Fz*uz)*
					(-0.2222222222222222 - (ux - uy)*(ux - uy) + 0.3333333333333333*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz))\
					- (m10 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((ux - uy)*(ux - uy)) + 0.3333333333333333*(2*ux + ux*ux - 2*uy + uy*uy + uz*uz)) + 
							0.0625*(Fx*(1 + ux) + Fy*(-1. + uy) + Fz*uz)*
							(-0.2222222222222222 - 1.*(ux*ux) + 2.*ux*uy - 1.*(uy*uy) + 
									0.3333333333333333*(2.*ux + ux*ux - 2.*uy + uy*uy + uz*uz)))/tau;

			// q = 11
			dist[11*Np+n] = m11 - 0.125*(Fx*(-1 + ux) + Fy*uy + Fz*(-1 + uz))*
					(0.2222222222222222 + (ux + uz)*(ux + uz) - 0.3333333333333333*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz))\
					- (m11 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((ux + uz)*(ux + uz)) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
							0.0625*(Fx*(-1. + ux) + Fy*uy + Fz*(-1. + uz))*
							(-0.2222222222222222 - 1.*(ux*ux) - 2.*ux*uz - 1.*(uz*uz) + 
									0.3333333333333333*(-2.*ux + ux*ux + uy*uy + (-2. + uz)*uz)))/tau;

			// q = 12
			dist[12*Np+n] = m12 + 0.125*(Fx + Fz + Fx*ux + Fy*uy + Fz*uz)*
					(-0.2222222222222222 - (ux + uz)*(ux + uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*(2 + uz)))\
					- (m12 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((ux + uz)*(ux + uz)) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
							0.0625*(Fx + Fz + Fx*ux + Fy*uy + Fz*uz)*
							(-0.2222222222222222 - 1.*(ux*ux) - 2.*ux*uz - 1.*(uz*uz) + 
									0.3333333333333333*(2.*ux + ux*ux + uy*uy + uz*(2. + uz))))/tau;

			// q = 13
			dist[13*Np+n] = m13 + 0.125*(Fz + Fx*(-1 + ux) + Fy*uy + Fz*uz)*
					(-0.2222222222222222 - (ux - uz)*(ux - uz) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*(2 + uz)))\
					- (m13 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((ux - uz)*(ux - uz)) + 0.3333333333333333*(-2*ux + ux*ux + uy*uy + uz*(2 + uz))) + 
							0.0625*(Fz + Fx*(-1. + ux) + Fy*uy + Fz*uz)*
							(-0.2222222222222222 - 1.*(ux*ux) + 2.*ux*uz - 1.*(uz*uz) + 
									0.3333333333333333*(-2.*ux + ux*ux + uy*uy + uz*(2. + uz))))/tau;

			// q= 14
			dist[14*Np+n] = m14 + 0.125*(Fx*(1 + ux) + Fy*uy + Fz*(-1 + uz))*
					(-0.2222222222222222 - (ux - uz)*(ux - uz) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz))\
					- (m14 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((ux - uz)*(ux - uz)) + 0.3333333333333333*(2*ux + ux*ux + uy*uy + (-2 + uz)*uz)) + 
							0.0625*(Fx*(1 + ux) + Fy*uy + Fz*(-1. + uz))*
							(-0.2222222222222222 - 1.*(ux*ux) + 2.*ux*uz - 1.*(uz*uz) + 
									0.3333333333333333*(2.*ux + ux*ux + uy*uy + (-2. + uz)*uz)))/tau;

			// q = 15
			dist[15*Np+n] = m15 - 0.125*(Fx*ux + Fy*(-1 + uy) + Fz*(-1 + uz))*
					(0.2222222222222222 + (uy + uz)*(uy + uz) - 0.3333333333333333*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz))\
					- (m15 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((uy + uz)*(uy + uz)) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + (-2 + uz)*uz)) + 
							0.0625*(Fx*ux + Fy*(-1. + uy) + Fz*(-1. + uz))*
							(-0.2222222222222222 - 1.*(uy*uy) - 2.*uy*uz - 1.*(uz*uz) + 
									0.3333333333333333*(ux*ux - 2.*uy + uy*uy + (-2. + uz)*uz)))/tau;

			// q = 16
			dist[16*Np+n] = m16 + 0.125*(Fy + Fz + Fx*ux + Fy*uy + Fz*uz)*
					(-0.2222222222222222 - (uy + uz)*(uy + uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*(2 + uz)))\
					- (m16 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((uy + uz)*(uy + uz)) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + uz*(2 + uz))) + 
							0.0625*(Fy + Fz + Fx*ux + Fy*uy + Fz*uz)*
							(-0.2222222222222222 - 1.*(uy*uy) - 2.*uy*uz - 1.*(uz*uz) + 
									0.3333333333333333*(ux*ux + 2.*uy + uy*uy + uz*(2. + uz))))/tau;

			// q = 17
			dist[17*Np+n] = m17 + 0.125*(Fz + Fx*ux + Fy*(-1 + uy) + Fz*uz)*
					(-0.2222222222222222 - (uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*(2 + uz)))\
					- (m17 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((uy - uz)*(uy - uz)) + 0.3333333333333333*(ux*ux - 2*uy + uy*uy + uz*(2 + uz))) + 
							0.0625*(Fz + Fx*ux + Fy*(-1. + uy) + Fz*uz)*
							(-0.2222222222222222 - 1.*(uy*uy) + 2.*uy*uz - 1.*(uz*uz) + 
									0.3333333333333333*(ux*ux - 2.*uy + uy*uy + uz*(2. + uz))))/tau;

			// q = 18
			dist[18*Np+n] = m18 + 0.125*(Fx*ux + Fy*(1 + uy) + Fz*(-1 + uz))*
					(-0.2222222222222222 - (uy - uz)*(uy - uz) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz))\
					- (m18 - 0.027777777777777776*p + 0.041666666666666664*rho0*
							(-((uy - uz)*(uy - uz)) + 0.3333333333333333*(ux*ux + 2*uy + uy*uy + (-2 + uz)*uz)) + 
							0.0625*(Fx*ux + Fy*(1 + uy) + Fz*(-1. + uz))*
							(-0.2222222222222222 - 1.*(uy*uy) + 2.*uy*uz - 1.*(uz*uz) + 
									0.3333333333333333*(ux*ux + 2.*uy + uy*uy + (-2. + uz)*uz)))/tau;
			//----------------------------------------------------------------------------------------------------------------------------------------//

			//Update velocity on device
			Vel[0*Np+n] = ux;
			Vel[1*Np+n] = uy;
			Vel[2*Np+n] = uz;
			//Update pressure on device
			Pressure[n] = p;
		}
	}
}

extern "C" void ScaLBL_D3Q19_FreeLeeModel_TwoFluid_Init(double *gqbar, double *mu_phi, double *ColorGrad, double Fx, double Fy, double Fz, int Np){

	dvc_ScaLBL_D3Q19_FreeLeeModel_TwoFluid_Init<<<NBLOCKS,NTHREADS >>>( gqbar,  mu_phi, ColorGrad, Fx, Fy, Fz, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_FreeLeeModel_TwoFluid_Init: %s \n",cudaGetErrorString(err));
	}
}


extern "C" void ScaLBL_D3Q19_FreeLeeModel_SingleFluid_Init(double *gqbar, double Fx, double Fy, double Fz, int Np){

	dvc_ScaLBL_D3Q19_FreeLeeModel_SingleFluid_Init<<<NBLOCKS,NTHREADS >>>( gqbar, Fx, Fy, Fz, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_FreeLeeModel_SingleFluid_Init: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_FreeLeeModel_PhaseField_Init(int *Map, double *Phi, double *Den, double *hq, double *ColorGrad, 
                                                    double rhoA, double rhoB, double tauM, double W, int start, int finish, int Np){
	
	dvc_ScaLBL_FreeLeeModel_PhaseField_Init<<<NBLOCKS,NTHREADS >>>(Map, Phi, Den, hq, ColorGrad, rhoA, rhoB, tauM, W, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_FreeLeeModel_PhaseField_Init: %s \n",cudaGetErrorString(err));
	}
	
	
}
extern "C" void ScaLBL_D3Q7_AAodd_FreeLee_PhaseField(int *neighborList, int *Map, double *hq, double *Den, double *Phi, double *ColorGrad, double *Vel,
                                                          double rhoA, double rhoB, double tauM, double W, int start, int finish, int Np)
{
//	cudaFuncSetCacheConfig(dvc_ScaLBL_D3Q7_AAodd_FreeLee_PhaseField, cudaFuncCachePreferL1);
//	dvc_ScaLBL_D3Q7_AAodd_FreeLee_PhaseField<<<NBLOCKS,NTHREADS >>>(neighborList, Map, hq, Den, Phi, ColorGrad, Vel,
 //            rhoA,  rhoB, tauM, W, start, finish,  Np);
//	cudaError_t err = cudaGetLastError();
//	if (cudaSuccess != err){
//		printf("CUDA error in ScaLBL_D3Q7_AAodd_FreeLee_PhaseField: %s \n",cudaGetErrorString(err));
//	}
}

extern "C" void ScaLBL_D3Q7_AAeven_FreeLee_PhaseField( int *Map, double *hq, double *Den, double *Phi, double *ColorGrad, double *Vel,
		double rhoA, double rhoB, double tauM, double W, int start, int finish, int Np){

//	cudaFuncSetCacheConfig(dvc_ScaLBL_D3Q7_AAeven_FreeLee_PhaseField, cudaFuncCachePreferL1);
//	dvc_ScaLBL_D3Q7_AAeven_FreeLee_PhaseField<<<NBLOCKS,NTHREADS >>>( Map, hq, Den, Phi, ColorGrad, Vel, rhoA, rhoB, tauM, W, start, finish, Np);
//	cudaError_t err = cudaGetLastError();
//	if (cudaSuccess != err){
//		printf("CUDA error in ScaLBL_D3Q7_AAeven_FreeLee_PhaseField: %s \n",cudaGetErrorString(err));
//	}
}

extern "C" void ScaLBL_D3Q7_AAodd_FreeLeeModel_PhaseField(int *neighborList, int *Map, double *hq, double *Den, double *Phi,
                                                          double rhoA, double rhoB, int start, int finish, int Np)
{
	cudaFuncSetCacheConfig(dvc_ScaLBL_D3Q7_AAodd_FreeLeeModel_PhaseField, cudaFuncCachePreferL1);
	dvc_ScaLBL_D3Q7_AAodd_FreeLeeModel_PhaseField<<<NBLOCKS,NTHREADS >>>(neighborList, Map, hq, Den, Phi, rhoA, rhoB, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_AAodd_FreeLeeModel_PhaseField: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q7_AAeven_FreeLeeModel_PhaseField( int *Map, double *hq, double *Den, double *Phi, 
		                                                    double rhoA, double rhoB, int start, int finish, int Np){

	cudaFuncSetCacheConfig(dvc_ScaLBL_D3Q7_AAeven_FreeLeeModel_PhaseField, cudaFuncCachePreferL1);
	dvc_ScaLBL_D3Q7_AAeven_FreeLeeModel_PhaseField<<<NBLOCKS,NTHREADS >>>( Map, hq, Den, Phi, rhoA, rhoB, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_AAeven_FreeLeeModel_PhaseField: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q7_ComputePhaseField(int *Map,  double *hq, double *Den, double *Phi, double rhoA, double rhoB, int start, int finish, int Np){

	cudaFuncSetCacheConfig(dvc_ScaLBL_D3Q7_ComputePhaseField, cudaFuncCachePreferL1);
	dvc_ScaLBL_D3Q7_ComputePhaseField<<<NBLOCKS,NTHREADS >>>( Map, hq, Den, Phi, rhoA, rhoB, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_ComputePhaseField: %s \n",cudaGetErrorString(err));
	}
}


extern "C" void ScaLBL_D3Q19_AAodd_FreeLeeModel(int *neighborList, int *Map, double *dist, double *Den,	double *Phi, double *mu_phi, double *Vel, double *Pressure, double *ColorGrad, 
                                                double rhoA, double rhoB, double tauA, double tauB, double kappa, double beta, double W, double Fx, double Fy, double Fz, 
                                                int strideY, int strideZ, int start, int finish, int Np){
	
	cudaFuncSetCacheConfig(dvc_ScaLBL_D3Q19_AAodd_FreeLeeModel, cudaFuncCachePreferL1);
	dvc_ScaLBL_D3Q19_AAodd_FreeLeeModel<<<NBLOCKS,NTHREADS >>>(neighborList, Map, dist, Den, Phi, mu_phi, Vel, Pressure,  ColorGrad, 
            rhoA, rhoB, tauA, tauB, kappa, beta, W, Fx, Fy, Fz, strideY, strideZ, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAodd_FreeLeeModel: %s \n",cudaGetErrorString(err));
	}
}	

extern "C" void ScaLBL_D3Q19_AAeven_FreeLeeModel(int *Map, double *dist, double *Den,	double *Phi, double *mu_phi, double *Vel, double *Pressure, double *ColorGrad,
                                                double rhoA, double rhoB, double tauA, double tauB, double kappa, double beta, double W, double Fx, double Fy, double Fz, 
                                                int strideY, int strideZ, int start, int finish, int Np){
	
	cudaFuncSetCacheConfig(dvc_ScaLBL_D3Q19_AAeven_FreeLeeModel, cudaFuncCachePreferL1);
	dvc_ScaLBL_D3Q19_AAeven_FreeLeeModel<<<NBLOCKS,NTHREADS >>>(Map, dist, Den, Phi, mu_phi, Vel, Pressure,  ColorGrad, 
            rhoA, rhoB, tauA, tauB, kappa, beta, W, Fx, Fy, Fz, strideY, strideZ, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_FreeLeeModel: %s \n",cudaGetErrorString(err));
	}
	
}

extern "C" void ScaLBL_D3Q19_AAodd_FreeLeeModel_Combined(int *neighborList, int *Map, double *dist, double *hq, double *Den,	double *Phi, double *mu_phi, double *Vel, double *Pressure, double *ColorGrad, 
                                                double rhoA, double rhoB, double tauA, double tauB, double tauM, double kappa, double beta, double W, double Fx, double Fy, double Fz, 
                                                int strideY, int strideZ, int start, int finish, int Np){
	
	cudaFuncSetCacheConfig(dvc_ScaLBL_D3Q19_AAodd_FreeLeeModel_Combined, cudaFuncCachePreferL1);
	dvc_ScaLBL_D3Q19_AAodd_FreeLeeModel_Combined<<<NBLOCKS,NTHREADS >>>(neighborList, Map, dist, hq, Den, Phi, mu_phi, Vel, Pressure,  ColorGrad, 
            rhoA, rhoB, tauA, tauB, tauM, kappa, beta, W, Fx, Fy, Fz, strideY, strideZ, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAodd_FreeLeeModel_Combined: %s \n",cudaGetErrorString(err));
	}
}	

extern "C" void ScaLBL_D3Q19_AAeven_FreeLeeModel_Combined(int *Map, double *dist, double *hq, double *Den,	double *Phi, double *mu_phi, double *Vel, double *Pressure, double *ColorGrad,
                                                double rhoA, double rhoB, double tauA, double tauB, double tauM, double kappa, double beta, double W, double Fx, double Fy, double Fz, 
                                                int strideY, int strideZ, int start, int finish, int Np){
	
	cudaFuncSetCacheConfig(dvc_ScaLBL_D3Q19_AAeven_FreeLeeModel_Combined, cudaFuncCachePreferL1);
	dvc_ScaLBL_D3Q19_AAeven_FreeLeeModel_Combined<<<NBLOCKS,NTHREADS >>>(Map, dist, hq, Den, Phi, mu_phi, Vel, Pressure,  ColorGrad, 
            rhoA, rhoB, tauA, tauB, tauM, kappa, beta, W, Fx, Fy, Fz, strideY, strideZ, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_FreeLeeModel_Combined: %s \n",cudaGetErrorString(err));
	}
	
}

extern "C" void ScaLBL_D3Q19_AAodd_FreeLeeModel_Combined_HigherOrder(int *neighborList, int *Map, double *dist, double *hq, double *Den,	double *Phi, double *mu_phi, double *Vel, double *Pressure, double *ColorGrad, 
                                                double rhoA, double rhoB, double tauA, double tauB, double tauM, double kappa, double beta, double W, double Fx, double Fy, double Fz, 
                                                int strideY, int strideZ, int start, int finish, int Np){
	
	cudaFuncSetCacheConfig(dvc_ScaLBL_D3Q19_AAodd_FreeLeeModel_Combined_HigherOrder, cudaFuncCachePreferL1);
	dvc_ScaLBL_D3Q19_AAodd_FreeLeeModel_Combined_HigherOrder<<<NBLOCKS,NTHREADS >>>(neighborList, Map, dist, hq, Den, Phi, mu_phi, Vel, Pressure,  ColorGrad, 
            rhoA, rhoB, tauA, tauB, tauM, kappa, beta, W, Fx, Fy, Fz, strideY, strideZ, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAodd_FreeLeeModel_Combined_HigherOrder: %s \n",cudaGetErrorString(err));
	}
}	

extern "C" void ScaLBL_D3Q19_AAeven_FreeLeeModel_Combined_HigherOrder(int *Map, double *dist, double *hq, double *Den,	double *Phi, double *mu_phi, double *Vel, double *Pressure, double *ColorGrad,
                                                double rhoA, double rhoB, double tauA, double tauB, double tauM, double kappa, double beta, double W, double Fx, double Fy, double Fz, 
                                                int strideY, int strideZ, int start, int finish, int Np){
	
	cudaFuncSetCacheConfig(dvc_ScaLBL_D3Q19_AAeven_FreeLeeModel_Combined_HigherOrder, cudaFuncCachePreferL1);
	dvc_ScaLBL_D3Q19_AAeven_FreeLeeModel_Combined_HigherOrder<<<NBLOCKS,NTHREADS >>>(Map, dist, hq, Den, Phi, mu_phi, Vel, Pressure,  ColorGrad, 
            rhoA, rhoB, tauA, tauB, tauM, kappa, beta, W, Fx, Fy, Fz, strideY, strideZ, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_FreeLeeModel_Combined_HigherOrder: %s \n",cudaGetErrorString(err));
	}
}
	
extern "C" void ScaLBL_D3Q19_AAodd_FreeLeeModel_SingleFluid_BGK(int *neighborList, double *dist, double *Vel, double *Pressure,  
                                                                double tau, double rho0, double Fx, double Fy, double Fz, int start, int finish, int Np){

	cudaFuncSetCacheConfig(dvc_ScaLBL_D3Q19_AAodd_FreeLeeModel_SingleFluid_BGK, cudaFuncCachePreferL1);
	dvc_ScaLBL_D3Q19_AAodd_FreeLeeModel_SingleFluid_BGK<<<NBLOCKS,NTHREADS >>>(neighborList, dist, Vel, Pressure, 
            tau, rho0, Fx, Fy, Fz, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAodd_FreeLeeModel_SingleFluid_BGK: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_AAeven_FreeLeeModel_SingleFluid_BGK(double *dist, double *Vel, double *Pressure, 
                                                                 double tau, double rho0, double Fx, double Fy, double Fz, int start, int finish, int Np){

	cudaFuncSetCacheConfig(dvc_ScaLBL_D3Q19_AAeven_FreeLeeModel_SingleFluid_BGK, cudaFuncCachePreferL1);
	dvc_ScaLBL_D3Q19_AAeven_FreeLeeModel_SingleFluid_BGK<<<NBLOCKS,NTHREADS >>>(dist, Vel, Pressure, 
            tau, rho0,  Fx, Fy, Fz, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_FreeLeeModel_SingleFluid_BGK: %s \n",cudaGetErrorString(err));
	}
}


extern "C" void ScaLBL_D3Q9_MGTest(int *Map, double *Phi,double *ColorGrad,int strideY, int strideZ, int start, int finish, int Np){
}
