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

#define NBLOCKS 1024
#define NTHREADS 256


#define CHECK_ERROR(KERNEL)                                         \
    do {                                                            \
        auto err = cudaGetLastError();                              \
        if ( cudaSuccess != err ){                                  \
            auto errString = cudaGetErrorString(err);               \
            printf("error in %s (kernel): %s \n",KERNEL,errString); \
        }                                                           \
    } while(0)


__global__ void dvc_ScaLBL_Solid_Dirichlet_D3Q7(double *dist, double *BoundaryValue, int *BounceBackDist_list, int *BounceBackSolid_list, int count)
{

    int idx;
    int iq,ib;
    double value_b,value_q;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		iq = BounceBackDist_list[idx];
        ib = BounceBackSolid_list[idx];
		value_b = BoundaryValue[ib];//get boundary value from a solid site
        value_q = dist[iq];
		dist[iq] = -1.0*value_q + value_b*0.25;//NOTE 0.25 is the speed of sound for D3Q7 lattice
	}
}

__global__ void dvc_ScaLBL_Solid_Neumann_D3Q7(double *dist, double *BoundaryValue, int *BounceBackDist_list, int *BounceBackSolid_list, int count)
{

    int idx;
    int iq,ib;
    double value_b,value_q;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		iq = BounceBackDist_list[idx];
        ib = BounceBackSolid_list[idx];
		value_b = BoundaryValue[ib];//get boundary value from a solid site
        value_q = dist[iq];
		dist[iq] = value_q + value_b;
	}
}

__global__ void dvc_ScaLBL_Solid_DirichletAndNeumann_D3Q7(double *dist, double *BoundaryValue,int *BoundaryLabel, int *BounceBackDist_list, int *BounceBackSolid_list, int count)
{

    int idx;
    int iq,ib;
    double value_b,value_b_label,value_q;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		iq = BounceBackDist_list[idx];
        ib = BounceBackSolid_list[idx];
		value_b = BoundaryValue[ib];//get boundary value from a solid site
		value_b_label = BoundaryLabel[ib];//get boundary label (i.e. type of BC) from a solid site
        value_q = dist[iq];
        if (value_b_label==1){//Dirichlet BC
		    dist[iq] = -1.0*value_q + value_b*0.25;//NOTE 0.25 is the speed of sound for D3Q7 lattice
        }
        if (value_b_label==2){//Neumann BC
		    dist[iq] = value_q + value_b;
        }
	}
}

__global__ void dvc_ScaLBL_Solid_SlippingVelocityBC_D3Q19(double *dist, double *zeta_potential, double *ElectricField, double *SolidGrad,
                                                          double epsilon_LB, double tau, double rho0,double den_scale, double h, double time_conv,
                                                          int *BounceBackDist_list, int *BounceBackSolid_list, int *FluidBoundary_list,
                                                          double *lattice_weight, float *lattice_cx, float *lattice_cy, float *lattice_cz,
                                                          int count, int Np)
{
    int idx;
    int iq,ib,ifluidBC;
    double value_b,value_q;
    double Ex,Ey,Ez;
    double Etx,Ety,Etz;//tangential part of electric field
    double E_mag_normal;
    double nsx,nsy,nsz;//unit normal solid gradient
    double ubx,uby,ubz;//slipping velocity at fluid boundary nodes
    float cx,cy,cz;//lattice velocity (D3Q19)
    double LB_weight;//lattice weighting coefficient (D3Q19)
    double cs2_inv = 3.0;//inverse of cs^2 for D3Q19
    double nu_LB = (tau-0.5)/cs2_inv;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		iq       = BounceBackDist_list[idx];
        ib       = BounceBackSolid_list[idx];
        ifluidBC = FluidBoundary_list[idx];
		value_b = zeta_potential[ib];//get zeta potential from a solid site
        value_q = dist[iq];

        //Load electric field and compute its tangential componet
        Ex = ElectricField[ifluidBC+0*Np]; 
        Ey = ElectricField[ifluidBC+1*Np];
        Ez = ElectricField[ifluidBC+2*Np];
        nsx = SolidGrad[ifluidBC+0*Np]; 
        nsy = SolidGrad[ifluidBC+1*Np];
        nsz = SolidGrad[ifluidBC+2*Np];
        E_mag_normal = Ex*nsx+Ey*nsy+Ez*nsz;//magnitude of electric field in the direction normal to solid nodes
        //compute tangential electric field
        Etx = Ex - E_mag_normal*nsx;
        Ety = Ey - E_mag_normal*nsy;
        Etz = Ez - E_mag_normal*nsz;
        ubx = -epsilon_LB*value_b*Etx/(nu_LB*rho0)*time_conv*time_conv/(h*h*1.0e-12)/den_scale;                                                                                                        
        uby = -epsilon_LB*value_b*Ety/(nu_LB*rho0)*time_conv*time_conv/(h*h*1.0e-12)/den_scale;                                                                                                        
        ubz = -epsilon_LB*value_b*Etz/(nu_LB*rho0)*time_conv*time_conv/(h*h*1.0e-12)/den_scale;                                                                                                        

        //compute bounce-back distribution
        LB_weight = lattice_weight[idx];
        cx = lattice_cx[idx];
        cy = lattice_cy[idx];
        cz = lattice_cz[idx];
		dist[iq] = value_q - 2.0*LB_weight*rho0*cs2_inv*(cx*ubx+cy*uby+cz*ubz);
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_z(int *list, double *dist, double Vin, int count, int Np)
{
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f6 = dist[5*Np+n];
		//...................................................
		f5 = Vin - (f0+f1+f2+f3+f4+f6);
		dist[6*Np+n] = f5;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_Z(int *list, double *dist, double Vout, int count, int Np)
{
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f5 = dist[6*Np+n];
		//...................................................
		f6 = Vout - (f0+f1+f2+f3+f4+f5);
		dist[5*Np+n] = f6;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_z(int *d_neighborList, int *list, double *dist, double Vin, int count, int Np)
{
	int idx, n;
    int nread,nr5;
	double f0,f1,f2,f3,f4,f5,f6;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];

		nread = d_neighborList[n];
		f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		f3 = dist[nread];

		nread = d_neighborList[n+Np];
		f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		f4 = dist[nread];

		nread = d_neighborList[n+5*Np];
		f6 = dist[nread];

		// Unknown distributions
		nr5 = d_neighborList[n+4*Np];
		f5 = Vin - (f0+f1+f2+f3+f4+f6);
		dist[nr5] = f5;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_Z(int *d_neighborList, int *list, double *dist, double Vout, int count, int Np)
{
	int idx, n;
    int nread,nr6;
	double f0,f1,f2,f3,f4,f5,f6;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];

		nread = d_neighborList[n];
		f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		f3 = dist[nread];

		nread = d_neighborList[n+4*Np];
		f5 = dist[nread];

		nread = d_neighborList[n+Np];
		f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		f4 = dist[nread];

		// unknown distributions
		nr6 = d_neighborList[n+5*Np];
		f6 = Vout - (f0+f1+f2+f3+f4+f5);
		dist[nr6] = f6;
	}
}

__global__ void dvc_ScaLBL_Poisson_D3Q7_BC_z(int *list, int *Map, double *Psi, double Vin, int count)
{
	int idx,n,nm;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		nm = Map[n];
		Psi[nm] = Vin;
	}
}


__global__ void dvc_ScaLBL_Poisson_D3Q7_BC_Z(int *list, int *Map, double *Psi, double Vout, int count)
{
	int idx,n,nm;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		nm = Map[n];
		Psi[nm] = Vout;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_z(int *list, double *dist, double Cin, int count, int Np)
{
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f6 = dist[5*Np+n];
		//...................................................
		f5 = Cin - (f0+f1+f2+f3+f4+f6);
		dist[6*Np+n] = f5;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_Z(int *list, double *dist, double Cout, int count, int Np)
{
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f5 = dist[6*Np+n];
		//...................................................
		f6 = Cout - (f0+f1+f2+f3+f4+f5);
		dist[5*Np+n] = f6;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_z(int *d_neighborList, int *list, double *dist, double Cin, int count, int Np)
{
	int idx, n;
    int nread,nr5;
	double f0,f1,f2,f3,f4,f5,f6;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];

		nread = d_neighborList[n];
		f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		f3 = dist[nread];

		nread = d_neighborList[n+Np];
		f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		f4 = dist[nread];

		nread = d_neighborList[n+5*Np];
		f6 = dist[nread];

		// Unknown distributions
		nr5 = d_neighborList[n+4*Np];
		f5 = Cin - (f0+f1+f2+f3+f4+f6);
		dist[nr5] = f5;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_Z(int *d_neighborList, int *list, double *dist, double Cout, int count, int Np)
{
	int idx, n;
    int nread,nr6;
	double f0,f1,f2,f3,f4,f5,f6;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];

		nread = d_neighborList[n];
		f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		f3 = dist[nread];

		nread = d_neighborList[n+4*Np];
		f5 = dist[nread];

		nread = d_neighborList[n+Np];
		f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		f4 = dist[nread];

		// unknown distributions
		nr6 = d_neighborList[n+5*Np];
		f6 = Cout - (f0+f1+f2+f3+f4+f5);
		dist[nr6] = f6;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_Diff_BC_z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np)
{
    //NOTE: FluxIn is the inward flux
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f6 = dist[5*Np+n];
        fsum_partial = f0+f1+f2+f3+f4+f6;
        uz = VelocityZ[n];
		//...................................................
        f5 =(FluxIn+(1.0-0.5/tau)*(f6+uz*fsum_partial))/(1.0-0.5/tau)/(1.0-uz); 
		dist[6*Np+n] = f5;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_Diff_BC_Z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np)
{
    //NOTE: FluxIn is the inward flux
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f5 = dist[6*Np+n];
        fsum_partial = f0+f1+f2+f3+f4+f5;
        uz = VelocityZ[n];
		//...................................................
        f6 =(FluxIn+(1.0-0.5/tau)*(f5-uz*fsum_partial))/(1.0-0.5/tau)/(1.0+uz); 
		dist[5*Np+n] = f6;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_Diff_BC_z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np)
{
    //NOTE: FluxIn is the inward flux
	int idx, n;
    int nread,nr5;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];

		nread = d_neighborList[n];
		f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		f3 = dist[nread];

		nread = d_neighborList[n+Np];
		f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		f4 = dist[nread];

		nread = d_neighborList[n+5*Np];
		f6 = dist[nread];

        fsum_partial = f0+f1+f2+f3+f4+f6;
        uz = VelocityZ[n];
		//...................................................
        f5 =(FluxIn+(1.0-0.5/tau)*(f6+uz*fsum_partial))/(1.0-0.5/tau)/(1.0-uz); 

		// Unknown distributions
		nr5 = d_neighborList[n+4*Np];
		dist[nr5] = f5;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_Diff_BC_Z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np)
{
    //NOTE: FluxIn is the inward flux
	int idx, n;
    int nread,nr6;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];

		nread = d_neighborList[n];
		f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		f3 = dist[nread];

		nread = d_neighborList[n+4*Np];
		f5 = dist[nread];

		nread = d_neighborList[n+Np];
		f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		f4 = dist[nread];

        fsum_partial = f0+f1+f2+f3+f4+f5;
        uz = VelocityZ[n];
		//...................................................
        f6 =(FluxIn+(1.0-0.5/tau)*(f5-uz*fsum_partial))/(1.0-0.5/tau)/(1.0+uz); 

		// unknown distributions
		nr6 = d_neighborList[n+5*Np];
		dist[nr6] = f6;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvc_BC_z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np)
{
    //NOTE: FluxIn is the inward flux
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f6 = dist[5*Np+n];
        fsum_partial = f0+f1+f2+f3+f4+f6;
        uz = VelocityZ[n];
		//...................................................
        f5 =(FluxIn+(1.0-0.5/tau)*f6-0.5*uz*fsum_partial/tau)/(1.0-0.5/tau+0.5*uz/tau); 
		dist[6*Np+n] = f5;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvc_BC_Z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np)
{
    //NOTE: FluxIn is the inward flux
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f5 = dist[6*Np+n];
        fsum_partial = f0+f1+f2+f3+f4+f5;
        uz = VelocityZ[n];
		//...................................................
        f6 =(FluxIn+(1.0-0.5/tau)*f5+0.5*uz*fsum_partial/tau)/(1.0-0.5/tau-0.5*uz/tau); 
		dist[5*Np+n] = f6;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvc_BC_z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np)
{
    //NOTE: FluxIn is the inward flux
	int idx, n;
    int nread,nr5;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];

		nread = d_neighborList[n];
		f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		f3 = dist[nread];

		nread = d_neighborList[n+Np];
		f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		f4 = dist[nread];

		nread = d_neighborList[n+5*Np];
		f6 = dist[nread];

        fsum_partial = f0+f1+f2+f3+f4+f6;
        uz = VelocityZ[n];
		//...................................................
        f5 =(FluxIn+(1.0-0.5/tau)*f6-0.5*uz*fsum_partial/tau)/(1.0-0.5/tau+0.5*uz/tau); 

		// Unknown distributions
		nr5 = d_neighborList[n+4*Np];
		dist[nr5] = f5;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvc_BC_Z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np)
{
    //NOTE: FluxIn is the inward flux
	int idx, n;
    int nread,nr6;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];

		nread = d_neighborList[n];
		f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		f3 = dist[nread];

		nread = d_neighborList[n+4*Np];
		f5 = dist[nread];

		nread = d_neighborList[n+Np];
		f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		f4 = dist[nread];

        fsum_partial = f0+f1+f2+f3+f4+f5;
        uz = VelocityZ[n];
		//...................................................
        f6 =(FluxIn+(1.0-0.5/tau)*f5+0.5*uz*fsum_partial/tau)/(1.0-0.5/tau-0.5*uz/tau); 

		// unknown distributions
		nr6 = d_neighborList[n+5*Np];
		dist[nr6] = f6;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvcElec_BC_z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, double *ElectricField_Z,
                                                                  double Di, double zi, double Vt, int count, int Np)
{
    //NOTE: FluxIn is the inward flux
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
    double uEPz;//electrochemical induced velocity
    double Ez;//electrical field
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f6 = dist[5*Np+n];
        fsum_partial = f0+f1+f2+f3+f4+f6;
        uz = VelocityZ[n];
        Ez = ElectricField_Z[n];
        uEPz=zi*Di/Vt*Ez;
		//...................................................
        f5 =(FluxIn+(1.0-0.5/tau)*f6-(0.5*uz/tau+uEPz)*fsum_partial)/(1.0-0.5/tau+0.5*uz/tau+uEPz); 
		dist[6*Np+n] = f5;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvcElec_BC_Z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, double *ElectricField_Z, 
                                                                  double Di, double zi, double Vt, int count, int Np)
{
    //NOTE: FluxIn is the inward flux
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
    double uEPz;//electrochemical induced velocity
    double Ez;//electrical field
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f5 = dist[6*Np+n];
        fsum_partial = f0+f1+f2+f3+f4+f5;
        uz = VelocityZ[n];
        Ez = ElectricField_Z[n];
        uEPz=zi*Di/Vt*Ez;
		//...................................................
        f6 =(FluxIn+(1.0-0.5/tau)*f5+(0.5*uz/tau+uEPz)*fsum_partial)/(1.0-0.5/tau-0.5*uz/tau-uEPz); 
		dist[5*Np+n] = f6;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvcElec_BC_z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, double *ElectricField_Z,
                                                                 double Di, double zi, double Vt, int count, int Np)
{
    //NOTE: FluxIn is the inward flux
	int idx, n;
    int nread,nr5;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
    double uEPz;//electrochemical induced velocity
    double Ez;//electrical field
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];

		nread = d_neighborList[n];
		f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		f3 = dist[nread];

		nread = d_neighborList[n+Np];
		f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		f4 = dist[nread];

		nread = d_neighborList[n+5*Np];
		f6 = dist[nread];

        fsum_partial = f0+f1+f2+f3+f4+f6;
        uz = VelocityZ[n];
        Ez = ElectricField_Z[n];
        uEPz=zi*Di/Vt*Ez;
		//...................................................
        f5 =(FluxIn+(1.0-0.5/tau)*f6-(0.5*uz/tau+uEPz)*fsum_partial)/(1.0-0.5/tau+0.5*uz/tau+uEPz); 

		// Unknown distributions
		nr5 = d_neighborList[n+4*Np];
		dist[nr5] = f5;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvcElec_BC_Z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, double *ElectricField_Z,
                                                                 double Di, double zi, double Vt, int count, int Np)
{
    //NOTE: FluxIn is the inward flux
	int idx, n;
    int nread,nr6;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
    double uEPz;//electrochemical induced velocity
    double Ez;//electrical field
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];

		nread = d_neighborList[n];
		f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		f3 = dist[nread];

		nread = d_neighborList[n+4*Np];
		f5 = dist[nread];

		nread = d_neighborList[n+Np];
		f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		f4 = dist[nread];

        fsum_partial = f0+f1+f2+f3+f4+f5;
        uz = VelocityZ[n];
        Ez = ElectricField_Z[n];
        uEPz=zi*Di/Vt*Ez;
		//...................................................
        f6 =(FluxIn+(1.0-0.5/tau)*f5+(0.5*uz/tau+uEPz)*fsum_partial)/(1.0-0.5/tau-0.5*uz/tau-uEPz); 

		// unknown distributions
		nr6 = d_neighborList[n+5*Np];
		dist[nr6] = f6;
	}
}
//*************************************************************************

extern "C" void ScaLBL_Solid_Dirichlet_D3Q7(double *dist, double *BoundaryValue, int *BounceBackDist_list, int *BounceBackSolid_list, int count){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_Solid_Dirichlet_D3Q7<<<GRID,512>>>(dist, BoundaryValue, BounceBackDist_list, BounceBackSolid_list, count);
    CHECK_ERROR("ScaLBL_Solid_Dirichlet_D3Q7");
}

extern "C" void ScaLBL_Solid_Neumann_D3Q7(double *dist, double *BoundaryValue, int *BounceBackDist_list, int *BounceBackSolid_list, int count){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_Solid_Neumann_D3Q7<<<GRID,512>>>(dist, BoundaryValue, BounceBackDist_list, BounceBackSolid_list, count);
    CHECK_ERROR("ScaLBL_Solid_Neumann_D3Q7");
}

extern "C" void ScaLBL_Solid_DirichletAndNeumann_D3Q7(double *dist, double *BoundaryValue,int *BoundaryLabel, int *BounceBackDist_list, int *BounceBackSolid_list, int count){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_Solid_DirichletAndNeumann_D3Q7<<<GRID,512>>>(dist, BoundaryValue, BoundaryLabel, BounceBackDist_list, BounceBackSolid_list, count);
    CHECK_ERROR("ScaLBL_Solid_DirichletAndNeumann_D3Q7");
}

extern "C" void ScaLBL_Solid_SlippingVelocityBC_D3Q19(double *dist, double *zeta_potential, double *ElectricField, double *SolidGrad,
                                                      double epsilon_LB, double tau, double rho0,double den_scale, double h, double time_conv,
                                                      int *BounceBackDist_list, int *BounceBackSolid_list, int *FluidBoundary_list,
                                                      double *lattice_weight, float *lattice_cx, float *lattice_cy, float *lattice_cz,
                                                      int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_Solid_SlippingVelocityBC_D3Q19<<<GRID,512>>>(dist, zeta_potential, ElectricField, SolidGrad,
                                                            epsilon_LB, tau, rho0, den_scale, h, time_conv,
                                                            BounceBackDist_list, BounceBackSolid_list, FluidBoundary_list,
                                                            lattice_weight, lattice_cx, lattice_cy, lattice_cz,
                                                            count, Np);
    CHECK_ERROR("ScaLBL_Solid_SlippingVelocityBC_D3Q19");
}

extern "C" void ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_z(int *list, double *dist, double Vin, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_z<<<GRID,512>>>(list, dist, Vin, count, Np);
    CHECK_ERROR("ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_z");
}

extern "C" void ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_Z(int *list, double *dist, double Vout, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_Z<<<GRID,512>>>(list, dist, Vout, count, Np);
    CHECK_ERROR("ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_Z");
}

extern "C" void ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_z(int *d_neighborList, int *list, double *dist, double Vin, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_z<<<GRID,512>>>(d_neighborList, list, dist, Vin, count, Np);
    CHECK_ERROR("ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_z");
}

extern "C" void ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_Z(int *d_neighborList, int *list, double *dist, double Vout, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_Z<<<GRID,512>>>(d_neighborList, list, dist, Vout, count, Np);
    CHECK_ERROR("ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_Z");
}

extern "C" void ScaLBL_Poisson_D3Q7_BC_z(int *list, int *Map, double *Psi, double Vin, int count){
	int GRID = count / 512 + 1;
    dvc_ScaLBL_Poisson_D3Q7_BC_z<<<GRID,512>>>(list, Map, Psi, Vin, count);
    CHECK_ERROR("ScaLBL_Poisson_D3Q7_BC_z");
}

extern "C" void ScaLBL_Poisson_D3Q7_BC_Z(int *list, int *Map, double *Psi, double Vout, int count){
	int GRID = count / 512 + 1;
    dvc_ScaLBL_Poisson_D3Q7_BC_Z<<<GRID,512>>>(list, Map, Psi, Vout, count);
    CHECK_ERROR("ScaLBL_Poisson_D3Q7_BC_Z");
}

extern "C" void ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_z(int *list, double *dist, double Cin, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_z<<<GRID,512>>>(list, dist, Cin, count, Np);
    CHECK_ERROR("ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_z");
}

extern "C" void ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_Z(int *list, double *dist, double Cout, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_Z<<<GRID,512>>>(list, dist, Cout, count, Np);
    CHECK_ERROR("ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_Z");
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_z(int *d_neighborList, int *list, double *dist, double Cin, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_z<<<GRID,512>>>(d_neighborList, list, dist, Cin, count, Np);
    CHECK_ERROR("ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_z");
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_Z(int *d_neighborList, int *list, double *dist, double Cout, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_Z<<<GRID,512>>>(d_neighborList, list, dist, Cout, count, Np);
    CHECK_ERROR("ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_Z");
}
//------------Diff-----------------
extern "C" void ScaLBL_D3Q7_AAeven_Ion_Flux_Diff_BC_z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_Diff_BC_z<<<GRID,512>>>(list, dist, FluxIn, tau, VelocityZ, count, Np);
    CHECK_ERROR("ScaLBL_D3Q7_AAeven_Ion_Flux_Diff_BC_z");
}

extern "C" void ScaLBL_D3Q7_AAeven_Ion_Flux_Diff_BC_Z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_Diff_BC_Z<<<GRID,512>>>(list, dist, FluxIn, tau, VelocityZ, count, Np);
    CHECK_ERROR("ScaLBL_D3Q7_AAeven_Ion_Flux_Diff_BC_Z");
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Flux_Diff_BC_z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_Diff_BC_z<<<GRID,512>>>(d_neighborList, list, dist, FluxIn, tau, VelocityZ, count, Np);
    CHECK_ERROR("ScaLBL_D3Q7_AAodd_Ion_Flux_Diff_BC_z");
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Flux_Diff_BC_Z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_Diff_BC_Z<<<GRID,512>>>(d_neighborList, list, dist, FluxIn, tau, VelocityZ, count, Np);
    CHECK_ERROR("ScaLBL_D3Q7_AAodd_Ion_Flux_Diff_BC_Z");
}
//----------DiffAdvc-------------
extern "C" void ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvc_BC_z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvc_BC_z<<<GRID,512>>>(list, dist, FluxIn, tau, VelocityZ, count, Np);
    CHECK_ERROR("ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvc_BC_z");
}

extern "C" void ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvc_BC_Z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvc_BC_Z<<<GRID,512>>>(list, dist, FluxIn, tau, VelocityZ, count, Np);
    CHECK_ERROR("ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvc_BC_Z");
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvc_BC_z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvc_BC_z<<<GRID,512>>>(d_neighborList, list, dist, FluxIn, tau, VelocityZ, count, Np);
    CHECK_ERROR("ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvc_BC_z");
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvc_BC_Z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvc_BC_Z<<<GRID,512>>>(d_neighborList, list, dist, FluxIn, tau, VelocityZ, count, Np);
    CHECK_ERROR("ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvc_BC_Z");
}
//----------DiffAdvcElec-------------
extern "C" void ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvcElec_BC_z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, double *ElectricField_Z,
                                                              double Di, double zi, double Vt, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvcElec_BC_z<<<GRID,512>>>(list, dist, FluxIn, tau, VelocityZ, ElectricField_Z, Di, zi, Vt, count, Np);
    CHECK_ERROR("ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvcElec_BC_z");
}

extern "C" void ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvcElec_BC_Z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, double *ElectricField_Z,
                                                              double Di, double zi, double Vt, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvcElec_BC_Z<<<GRID,512>>>(list, dist, FluxIn, tau, VelocityZ, ElectricField_Z, Di, zi, Vt, count, Np);
    CHECK_ERROR("ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvcElec_BC_Z");
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvcElec_BC_z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, double *ElectricField_Z,
                                                             double Di, double zi, double Vt, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvcElec_BC_z<<<GRID,512>>>(d_neighborList, list, dist, FluxIn, tau, VelocityZ, ElectricField_Z, Di, zi, Vt, count, Np);
    CHECK_ERROR("ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvcElec_BC_z");
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvcElec_BC_Z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, double *ElectricField_Z,
                                                             double Di, double zi, double Vt, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvcElec_BC_Z<<<GRID,512>>>(d_neighborList, list, dist, FluxIn, tau, VelocityZ, ElectricField_Z, Di, zi, Vt, count, Np);
    CHECK_ERROR("ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvcElec_BC_Z");
}
//-------------------------------
