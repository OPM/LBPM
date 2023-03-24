#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>
#include <math.h>
#include <stdio.h>

#define NBLOCKS 1024
#define NTHREADS 256

/*
DPCT1010:218: SYCL uses exceptions to report errors and does not use the error
codes. The call was replaced with 0. You need to rewrite this code.
*/
/*
DPCT1009:219: SYCL uses exceptions to report errors and does not use the error
codes. The original code was commented out and a warning string was inserted.
You need to rewrite this code.
*/
#define CHECK_ERROR(KERNEL)                                                    \
    do {                                                                       \
        auto err = 0;                                                          \
                                                                               \
    } while (0)

void dvc_ScaLBL_Solid_Dirichlet_D3Q7(double *dist, double *BoundaryValue, int *BounceBackDist_list, int *BounceBackSolid_list, int count,
                                     sycl::nd_item<3> item_ct1)
{

    int idx;
    int iq,ib;
    double value_b,value_q;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
        if (idx < count){
		iq = BounceBackDist_list[idx];
        ib = BounceBackSolid_list[idx];
		value_b = BoundaryValue[ib];//get boundary value from a solid site
        value_q = dist[iq];
		dist[iq] = -1.0*value_q + value_b*0.25;//NOTE 0.25 is the speed of sound for D3Q7 lattice
	}
}

void dvc_ScaLBL_Solid_Neumann_D3Q7(double *dist, double *BoundaryValue, int *BounceBackDist_list, int *BounceBackSolid_list, int count,
                                   sycl::nd_item<3> item_ct1)
{

    int idx;
    int iq,ib;
    double value_b,value_q;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
        if (idx < count){
		iq = BounceBackDist_list[idx];
        ib = BounceBackSolid_list[idx];
		value_b = BoundaryValue[ib];//get boundary value from a solid site
        value_q = dist[iq];
		dist[iq] = value_q + value_b;
	}
}

void dvc_ScaLBL_Solid_DirichletAndNeumann_D3Q7(double *dist, double *BoundaryValue,int *BoundaryLabel, int *BounceBackDist_list, int *BounceBackSolid_list, int count,
                                               sycl::nd_item<3> item_ct1)
{

    int idx;
    int iq,ib;
    double value_b,value_b_label,value_q;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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

void dvc_ScaLBL_Solid_SlippingVelocityBC_D3Q19(double *dist, double *zeta_potential, double *ElectricField, double *SolidGrad,
                                                          double epsilon_LB, double tau, double rho0,double den_scale, double h, double time_conv,
                                                          int *BounceBackDist_list, int *BounceBackSolid_list, int *FluidBoundary_list,
                                                          double *lattice_weight, float *lattice_cx, float *lattice_cy, float *lattice_cz,
                                                          int count, int Np,
                                                          sycl::nd_item<3> item_ct1)
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
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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

void dvc_ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_z(int *list, double *dist, double Vin, int count, int Np,
                                                   sycl::nd_item<3> item_ct1)
{
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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

void dvc_ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_Z(int *list, double *dist, double Vout, int count, int Np,
                                                   sycl::nd_item<3> item_ct1)
{
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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

void dvc_ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_z(int *d_neighborList, int *list, double *dist, double Vin, int count, int Np,
                                                  sycl::nd_item<3> item_ct1)
{
	int idx, n;
    int nread,nr5;
	double f0,f1,f2,f3,f4,f5,f6;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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

void dvc_ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_Z(int *d_neighborList, int *list, double *dist, double Vout, int count, int Np,
                                                  sycl::nd_item<3> item_ct1)
{
	int idx, n;
    int nread,nr6;
	double f0,f1,f2,f3,f4,f5,f6;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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

void dvc_ScaLBL_Poisson_D3Q7_BC_z(int *list, int *Map, double *Psi, double Vin, int count,
                                  sycl::nd_item<3> item_ct1)
{
	int idx,n,nm;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
        if (idx < count){
		n = list[idx];
		nm = Map[n];
		Psi[nm] = Vin;
	}
}


void dvc_ScaLBL_Poisson_D3Q7_BC_Z(int *list, int *Map, double *Psi, double Vout, int count,
                                  sycl::nd_item<3> item_ct1)
{
	int idx,n,nm;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
        if (idx < count){
		n = list[idx];
		nm = Map[n];
		Psi[nm] = Vout;
	}
}

void dvc_ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_z(int *list, double *dist, double Cin, int count, int Np,
                                                   sycl::nd_item<3> item_ct1)
{
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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

void dvc_ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_Z(int *list, double *dist, double Cout, int count, int Np,
                                                   sycl::nd_item<3> item_ct1)
{
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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

void dvc_ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_z(int *d_neighborList, int *list, double *dist, double Cin, int count, int Np,
                                                  sycl::nd_item<3> item_ct1)
{
	int idx, n;
    int nread,nr5;
	double f0,f1,f2,f3,f4,f5,f6;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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

void dvc_ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_Z(int *d_neighborList, int *list, double *dist, double Cout, int count, int Np,
                                                  sycl::nd_item<3> item_ct1)
{
	int idx, n;
    int nread,nr6;
	double f0,f1,f2,f3,f4,f5,f6;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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

void dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_Diff_BC_z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np,
                                               sycl::nd_item<3> item_ct1)
{
    //NOTE: FluxIn is the inward flux
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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

void dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_Diff_BC_Z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np,
                                               sycl::nd_item<3> item_ct1)
{
    //NOTE: FluxIn is the inward flux
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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

void dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_Diff_BC_z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np,
                                              sycl::nd_item<3> item_ct1)
{
    //NOTE: FluxIn is the inward flux
	int idx, n;
    int nread,nr5;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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

void dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_Diff_BC_Z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np,
                                              sycl::nd_item<3> item_ct1)
{
    //NOTE: FluxIn is the inward flux
	int idx, n;
    int nread,nr6;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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

void dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvc_BC_z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np,
                                                   sycl::nd_item<3> item_ct1)
{
    //NOTE: FluxIn is the inward flux
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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

void dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvc_BC_Z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np,
                                                   sycl::nd_item<3> item_ct1)
{
    //NOTE: FluxIn is the inward flux
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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

void dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvc_BC_z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np,
                                                  sycl::nd_item<3> item_ct1)
{
    //NOTE: FluxIn is the inward flux
	int idx, n;
    int nread,nr5;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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

void dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvc_BC_Z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np,
                                                  sycl::nd_item<3> item_ct1)
{
    //NOTE: FluxIn is the inward flux
	int idx, n;
    int nread,nr6;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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

void dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvcElec_BC_z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, double *ElectricField_Z,
                                                                  double Di, double zi, double Vt, int count, int Np,
                                                                  sycl::nd_item<3> item_ct1)
{
    //NOTE: FluxIn is the inward flux
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
    double uEPz;//electrochemical induced velocity
    double Ez;//electrical field
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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

void dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvcElec_BC_Z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, double *ElectricField_Z, 
                                                                  double Di, double zi, double Vt, int count, int Np,
                                                                  sycl::nd_item<3> item_ct1)
{
    //NOTE: FluxIn is the inward flux
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
    double uEPz;//electrochemical induced velocity
    double Ez;//electrical field
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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

void dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvcElec_BC_z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, double *ElectricField_Z,
                                                                 double Di, double zi, double Vt, int count, int Np,
                                                                 sycl::nd_item<3> item_ct1)
{
    //NOTE: FluxIn is the inward flux
	int idx, n;
    int nread,nr5;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
    double uEPz;//electrochemical induced velocity
    double Ez;//electrical field
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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

void dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvcElec_BC_Z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, double *ElectricField_Z,
                                                                 double Di, double zi, double Vt, int count, int Np,
                                                                 sycl::nd_item<3> item_ct1)
{
    //NOTE: FluxIn is the inward flux
	int idx, n;
    int nread,nr6;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
    double uEPz;//electrochemical induced velocity
    double Ez;//electrical field
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
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
        /*
        DPCT1049:45: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_Solid_Dirichlet_D3Q7(
                        dist, BoundaryValue, BounceBackDist_list,
                        BounceBackSolid_list, count, item_ct1);
            });
    CHECK_ERROR("ScaLBL_Solid_Dirichlet_D3Q7");
}

extern "C" void ScaLBL_Solid_Neumann_D3Q7(double *dist, double *BoundaryValue, int *BounceBackDist_list, int *BounceBackSolid_list, int count){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:46: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_Solid_Neumann_D3Q7(
                        dist, BoundaryValue, BounceBackDist_list,
                        BounceBackSolid_list, count, item_ct1);
            });
    CHECK_ERROR("ScaLBL_Solid_Neumann_D3Q7");
}

extern "C" void ScaLBL_Solid_DirichletAndNeumann_D3Q7(double *dist, double *BoundaryValue,int *BoundaryLabel, int *BounceBackDist_list, int *BounceBackSolid_list, int count){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:47: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_Solid_DirichletAndNeumann_D3Q7(
                        dist, BoundaryValue, BoundaryLabel, BounceBackDist_list,
                        BounceBackSolid_list, count, item_ct1);
            });
    CHECK_ERROR("ScaLBL_Solid_DirichletAndNeumann_D3Q7");
}

extern "C" void ScaLBL_Solid_SlippingVelocityBC_D3Q19(double *dist, double *zeta_potential, double *ElectricField, double *SolidGrad,
                                                      double epsilon_LB, double tau, double rho0,double den_scale, double h, double time_conv,
                                                      int *BounceBackDist_list, int *BounceBackSolid_list, int *FluidBoundary_list,
                                                      double *lattice_weight, float *lattice_cx, float *lattice_cy, float *lattice_cz,
                                                      int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:48: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_Solid_SlippingVelocityBC_D3Q19(
                        dist, zeta_potential, ElectricField, SolidGrad,
                        epsilon_LB, tau, rho0, den_scale, h, time_conv,
                        BounceBackDist_list, BounceBackSolid_list,
                        FluidBoundary_list, lattice_weight, lattice_cx,
                        lattice_cy, lattice_cz, count, Np, item_ct1);
            });
    CHECK_ERROR("ScaLBL_Solid_SlippingVelocityBC_D3Q19");
}

extern "C" void ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_z(int *list, double *dist, double Vin, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:49: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_z(
                        list, dist, Vin, count, Np, item_ct1);
            });
    CHECK_ERROR("ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_z");
}

extern "C" void ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_Z(int *list, double *dist, double Vout, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:50: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_Z(
                        list, dist, Vout, count, Np, item_ct1);
            });
    CHECK_ERROR("ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_Z");
}

extern "C" void ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_z(int *d_neighborList, int *list, double *dist, double Vin, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:51: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_z(
                        d_neighborList, list, dist, Vin, count, Np, item_ct1);
            });
    CHECK_ERROR("ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_z");
}

extern "C" void ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_Z(int *d_neighborList, int *list, double *dist, double Vout, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:52: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_Z(
                        d_neighborList, list, dist, Vout, count, Np, item_ct1);
            });
    CHECK_ERROR("ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_Z");
}

extern "C" void ScaLBL_Poisson_D3Q7_BC_z(int *list, int *Map, double *Psi, double Vin, int count){
	int GRID = count / 512 + 1;
    /*
    DPCT1049:53: The work-group size passed to the SYCL kernel may exceed the
    limit. To get the device limit, query info::device::max_work_group_size.
    Adjust the work-group size if needed.
    */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_Poisson_D3Q7_BC_z(list, Map, Psi, Vin, count,
                                                 item_ct1);
            });
    CHECK_ERROR("ScaLBL_Poisson_D3Q7_BC_z");
}

extern "C" void ScaLBL_Poisson_D3Q7_BC_Z(int *list, int *Map, double *Psi, double Vout, int count){
	int GRID = count / 512 + 1;
    /*
    DPCT1049:54: The work-group size passed to the SYCL kernel may exceed the
    limit. To get the device limit, query info::device::max_work_group_size.
    Adjust the work-group size if needed.
    */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_Poisson_D3Q7_BC_Z(list, Map, Psi, Vout, count,
                                                 item_ct1);
            });
    CHECK_ERROR("ScaLBL_Poisson_D3Q7_BC_Z");
}

extern "C" void ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_z(int *list, double *dist, double Cin, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:55: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_z(
                        list, dist, Cin, count, Np, item_ct1);
            });
    CHECK_ERROR("ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_z");
}

extern "C" void ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_Z(int *list, double *dist, double Cout, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:56: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_Z(
                        list, dist, Cout, count, Np, item_ct1);
            });
    CHECK_ERROR("ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_Z");
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_z(int *d_neighborList, int *list, double *dist, double Cin, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:57: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_z(
                        d_neighborList, list, dist, Cin, count, Np, item_ct1);
            });
    CHECK_ERROR("ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_z");
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_Z(int *d_neighborList, int *list, double *dist, double Cout, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:58: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_Z(
                        d_neighborList, list, dist, Cout, count, Np, item_ct1);
            });
    CHECK_ERROR("ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_Z");
}
//------------Diff-----------------
extern "C" void ScaLBL_D3Q7_AAeven_Ion_Flux_Diff_BC_z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:59: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_Diff_BC_z(
                        list, dist, FluxIn, tau, VelocityZ, count, Np,
                        item_ct1);
            });
    CHECK_ERROR("ScaLBL_D3Q7_AAeven_Ion_Flux_Diff_BC_z");
}

extern "C" void ScaLBL_D3Q7_AAeven_Ion_Flux_Diff_BC_Z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:60: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_Diff_BC_Z(
                        list, dist, FluxIn, tau, VelocityZ, count, Np,
                        item_ct1);
            });
    CHECK_ERROR("ScaLBL_D3Q7_AAeven_Ion_Flux_Diff_BC_Z");
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Flux_Diff_BC_z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:61: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_Diff_BC_z(
                        d_neighborList, list, dist, FluxIn, tau, VelocityZ,
                        count, Np, item_ct1);
            });
    CHECK_ERROR("ScaLBL_D3Q7_AAodd_Ion_Flux_Diff_BC_z");
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Flux_Diff_BC_Z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:62: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_Diff_BC_Z(
                        d_neighborList, list, dist, FluxIn, tau, VelocityZ,
                        count, Np, item_ct1);
            });
    CHECK_ERROR("ScaLBL_D3Q7_AAodd_Ion_Flux_Diff_BC_Z");
}
//----------DiffAdvc-------------
extern "C" void ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvc_BC_z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:63: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvc_BC_z(
                        list, dist, FluxIn, tau, VelocityZ, count, Np,
                        item_ct1);
            });
    CHECK_ERROR("ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvc_BC_z");
}

extern "C" void ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvc_BC_Z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:64: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvc_BC_Z(
                        list, dist, FluxIn, tau, VelocityZ, count, Np,
                        item_ct1);
            });
    CHECK_ERROR("ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvc_BC_Z");
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvc_BC_z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:65: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvc_BC_z(
                        d_neighborList, list, dist, FluxIn, tau, VelocityZ,
                        count, Np, item_ct1);
            });
    CHECK_ERROR("ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvc_BC_z");
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvc_BC_Z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:66: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvc_BC_Z(
                        d_neighborList, list, dist, FluxIn, tau, VelocityZ,
                        count, Np, item_ct1);
            });
    CHECK_ERROR("ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvc_BC_Z");
}
//----------DiffAdvcElec-------------
extern "C" void ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvcElec_BC_z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, double *ElectricField_Z,
                                                              double Di, double zi, double Vt, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:67: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvcElec_BC_z(
                        list, dist, FluxIn, tau, VelocityZ, ElectricField_Z, Di,
                        zi, Vt, count, Np, item_ct1);
            });
    CHECK_ERROR("ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvcElec_BC_z");
}

extern "C" void ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvcElec_BC_Z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, double *ElectricField_Z,
                                                              double Di, double zi, double Vt, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:68: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvcElec_BC_Z(
                        list, dist, FluxIn, tau, VelocityZ, ElectricField_Z, Di,
                        zi, Vt, count, Np, item_ct1);
            });
    CHECK_ERROR("ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvcElec_BC_Z");
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvcElec_BC_z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, double *ElectricField_Z,
                                                             double Di, double zi, double Vt, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:69: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvcElec_BC_z(
                        d_neighborList, list, dist, FluxIn, tau, VelocityZ,
                        ElectricField_Z, Di, zi, Vt, count, Np, item_ct1);
            });
    CHECK_ERROR("ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvcElec_BC_z");
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvcElec_BC_Z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, double *ElectricField_Z,
                                                             double Di, double zi, double Vt, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:70: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvcElec_BC_Z(
                        d_neighborList, list, dist, FluxIn, tau, VelocityZ,
                        ElectricField_Z, Di, zi, Vt, count, Np, item_ct1);
            });
    CHECK_ERROR("ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvcElec_BC_Z");
}
//-------------------------------
