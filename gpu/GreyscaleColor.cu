#include <stdio.h>

#define NBLOCKS 1024
#define NTHREADS 256

__global__ void dvc_ScaLBL_D3Q19_GreyscaleColor_Pressure(double *dist, double *Den, double *Poros,double *Velocity,
                double *Pressure, double rhoA,double rhoB, int N){

	int n;
    double ux,uy,uz,u_mag;
    double pressure;
    double porosity;
    double rho0;
    double phi;
    double nA,nB;

	int S = N/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
	    //........Get 1-D index for this thread....................
	    n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;

		if (n<N){		

            // initialize pressure value
            pressure = 0.0;
            pressure +=dist[1*N+n];
            pressure +=dist[2*N+n];
            pressure +=dist[3*N+n];
            pressure +=dist[4*N+n];
            pressure +=dist[5*N+n];
            pressure +=dist[6*N+n];
            pressure +=dist[7*N+n];
            pressure +=dist[8*N+n];
            pressure +=dist[9*N+n];
            pressure +=dist[10*N+n];
            pressure +=dist[11*N+n];
            pressure +=dist[12*N+n];
            pressure +=dist[13*N+n];
            pressure +=dist[14*N+n];
            pressure +=dist[15*N+n];
            pressure +=dist[16*N+n];
            pressure +=dist[17*N+n];
            pressure +=dist[18*N+n];

			// read the component number densities
			nA = Den[n];
			nB = Den[N + n];
			// compute phase indicator field
			phi=(nA-nB)/(nA+nB);
			// local density
			rho0=rhoA + 0.5*(1.0-phi)*(rhoB-rhoA);
            // read voxel porosity 
            porosity = Poros[n];
            // read velocity
            ux = Velocity[0*N+n]; 
            uy = Velocity[1*N+n];
            uz = Velocity[2*N+n];
            u_mag=sqrt(ux*ux+uy*uy+uz*uz);

            //Calculate pressure for Incompressible-MRT model
            pressure=0.5/porosity*(pressure-0.5*rho0*u_mag*u_mag/porosity);
            //Update pressure on device
            Pressure[n] = pressure;
		}
	}
}


__global__ void dvc_ScaLBL_D3Q19_AAeven_GreyscaleColor(double *dist, double *Aq, double *Bq, double *Den,
                double *DenGradA, double *DenGradB, double *SolidForce, int start, int finish, int Np,
                double tauA,double tauB,double tauA_eff,double tauB_eff,double rhoA,double rhoB,double Gsc, double Gx, double Gy, double Gz,
                double *Poros,double *Perm, double *Velocity,double *Pressure){

	int n;
	double vx,vy,vz,v_mag;
    double ux,uy,uz,u_mag;
    double pressure;//defined for this incompressible model
	// conserved momemnts
	double jx,jy,jz;
	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
    double fq;
    // currently disable 'GeoFun'
    double GeoFun=0.0;//geometric function from Guo's PRE 66, 036304 (2002)
    double porosity;
    double perm;//voxel permeability
    double c0, c1; //Guo's model parameters
    double Fx, Fy, Fz;//The total body force including Brinkman force and user-specified (Gx,Gy,Gz)
	double tau,tau_eff,rlx_setA,rlx_setB;
    double mu_eff;//effective kinematic viscosity for Darcy term
    double rho0;
    double phi;
    double nx,ny,nz,C;
    double nA,nB;
	double a1,b1,a2,b2,nAB,delta;
    double beta=0.95;
    double nA_gradx,nA_grady,nA_gradz;
    double nB_gradx,nB_grady,nB_gradz;
    double Gff_x,Gff_y,Gff_z;
    double Gfs_x,Gfs_y,Gfs_z;


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


	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
	    //........Get 1-D index for this thread....................
	    n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

		if ( n<finish ){

			// read the component number densities
			nA = Den[n];
			nB = Den[Np + n];
			nA_gradx = DenGradA[n+0*Np];
			nA_grady = DenGradA[n+1*Np];
			nA_gradz = DenGradA[n+2*Np];
			nB_gradx = DenGradB[n+0*Np];
			nB_grady = DenGradB[n+1*Np];
			nB_gradz = DenGradB[n+2*Np];

			// compute phase indicator field
			phi=(nA-nB)/(nA+nB);

			// local density
			rho0=rhoA + 0.5*(1.0-phi)*(rhoB-rhoA);
			// local relaxation time
			tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);
			rlx_setA = 1.f/tau;
			rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
			tau_eff=tauA_eff + 0.5*(1.0-phi)*(tauB_eff-tauA_eff);
            mu_eff = (tau_eff-0.5)/3.f;//kinematic viscosity


            //........................................................................
            //					READ THE DISTRIBUTIONS
            //		(read from opposite array due to previous swap operation)
            //........................................................................
            // q=0
            fq = dist[n];
            m1  = -30.0*fq;
            m2  = 12.0*fq;

            // q=1
            fq = dist[2*Np+n];
            pressure = fq;
            m1 -= 11.0*fq;
            m2 -= 4.0*fq;
            jx = fq;
            m4 = -4.0*fq;
            m9 = 2.0*fq;
            m10 = -4.0*fq;

            // f2 = dist[10*Np+n];
            fq = dist[1*Np+n];
            pressure += fq;
            m1 -= 11.0*(fq);
            m2 -= 4.0*(fq);
            jx -= fq;
            m4 += 4.0*(fq);
            m9 += 2.0*(fq);
            m10 -= 4.0*(fq);

            // q=3
            fq = dist[4*Np+n];
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            //---------------------------------------------------------------------//

            //---------------- Calculate SC fluid-fluid and fluid-solid forces ---------------//
            // fluid-fluid force
//            Gff_x = -Gsc*nA*nB_gradx*int(phi>0.0)-Gsc*nB*nA_gradx*int(phi<0.0);
//            Gff_y = -Gsc*nA*nB_grady*int(phi>0.0)-Gsc*nB*nA_grady*int(phi<0.0);
//            Gff_z = -Gsc*nA*nB_gradz*int(phi>0.0)-Gsc*nB*nA_gradz*int(phi<0.0);
            Gff_x = -Gsc*(nA*nB_gradx+nB*nA_gradx);
            Gff_y = -Gsc*(nA*nB_grady+nB*nA_grady);
            Gff_z = -Gsc*(nA*nB_gradz+nB*nA_gradz);
            // fluid-solid force
            Gfs_x = (nA-nB)*SolidForce[n+0*Np];    
            Gfs_y = (nA-nB)*SolidForce[n+1*Np];    
            Gfs_z = (nA-nB)*SolidForce[n+2*Np];    

            porosity = Poros[n];
            // use local saturation as an estimation of effective relperm values
            perm = Perm[n]*nA/(nA+nB)*int(phi>0.0)+Perm[n]*nB/(nA+nB)*int(phi<0.0);

            c0 = 0.5*(1.0+porosity*0.5*mu_eff/perm);
            if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
            //GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
            c1 = porosity*0.5*GeoFun/sqrt(perm);
            if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes

            vx = jx/rho0+0.5*(porosity*Gx+Gff_x+Gfs_x);
            vy = jy/rho0+0.5*(porosity*Gy+Gff_y+Gfs_y);
            vz = jz/rho0+0.5*(porosity*Gz+Gff_z+Gfs_z);
            v_mag=sqrt(vx*vx+vy*vy+vz*vz);
            ux = vx/(c0+sqrt(c0*c0+c1*v_mag));
            uy = vy/(c0+sqrt(c0*c0+c1*v_mag));
            uz = vz/(c0+sqrt(c0*c0+c1*v_mag));
            u_mag=sqrt(ux*ux+uy*uy+uz*uz);

            //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
            Fx = rho0*(-porosity*mu_eff/perm*ux - porosity*GeoFun/sqrt(perm)*u_mag*ux + porosity*Gx + Gff_x + Gfs_x);
            Fy = rho0*(-porosity*mu_eff/perm*uy - porosity*GeoFun/sqrt(perm)*u_mag*uy + porosity*Gy + Gff_y + Gfs_y);
            Fz = rho0*(-porosity*mu_eff/perm*uz - porosity*GeoFun/sqrt(perm)*u_mag*uz + porosity*Gz + Gff_z + Gfs_z);
            if (porosity==1.0){
                Fx=rho0*(Gx + Gff_x + Gfs_x);
                Fy=rho0*(Gy + Gff_y + Gfs_y);
                Fz=rho0*(Gz + Gff_z + Gfs_z);
            }

            //Calculate pressure for Incompressible-MRT model
            pressure=0.5/porosity*(pressure-0.5*rho0*u_mag*u_mag/porosity);

//            //..............carry out relaxation process...............................................
//            m1 = m1 + rlx_setA*((-30*rho0+19*(ux*ux+uy*uy+uz*uz)/porosity + 57*pressure*porosity) - m1) 
//                    + (1-0.5*rlx_setA)*38*(Fx*ux+Fy*uy+Fz*uz)/porosity;
//            m2 = m2 + rlx_setA*((12*rho0 - 5.5*(ux*ux+uy*uy+uz*uz)/porosity-27*pressure*porosity) - m2)
//                    + (1-0.5*rlx_setA)*11*(-Fx*ux-Fy*uy-Fz*uz)/porosity;
//            jx = jx + Fx;
//            m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0) - m4)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
//            jy = jy + Fy;
//            m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0) - m6)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
//            jz = jz + Fz;
//            m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0) - m8)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
//            m9 = m9 + rlx_setA*((rho0*(2*ux*ux-uy*uy-uz*uz)/porosity) - m9)
//                    + (1-0.5*rlx_setA)*(4*Fx*ux-2*Fy*uy-2*Fz*uz)/porosity;
//            m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10)
//                      + (1-0.5*rlx_setA)*(-2*Fx*ux+Fy*uy+Fz*uz)/porosity;
//            m11 = m11 + rlx_setA*((rho0*(uy*uy-uz*uz)/porosity) - m11)
//                      + (1-0.5*rlx_setA)*(2*Fy*uy-2*Fz*uz)/porosity;
//            m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12)
//                      + (1-0.5*rlx_setA)*(-Fy*uy+Fz*uz)/porosity;
//            m13 = m13 + rlx_setA*((rho0*ux*uy/porosity) - m13)
//                      + (1-0.5*rlx_setA)*(Fy*ux+Fx*uy)/porosity;
//            m14 = m14 + rlx_setA*((rho0*uy*uz/porosity) - m14)
//                      + (1-0.5*rlx_setA)*(Fz*uy+Fy*uz)/porosity;
//            m15 = m15 + rlx_setA*((rho0*ux*uz/porosity) - m15)
//                      + (1-0.5*rlx_setA)*(Fz*ux+Fx*uz)/porosity;
//            m16 = m16 + rlx_setB*( - m16);
//            m17 = m17 + rlx_setB*( - m17);
//            m18 = m18 + rlx_setB*( - m18);
//            //.......................................................................................................

            //-------------------- IMRT collison where body force has NO higher-order terms -------------//
            //..............carry out relaxation process...............................................
            m1 = m1 + rlx_setA*((-30*rho0+19*(ux*ux+uy*uy+uz*uz)/porosity + 57*pressure*porosity) - m1);
            m2 = m2 + rlx_setA*((12*rho0 - 5.5*(ux*ux+uy*uy+uz*uz)/porosity-27*pressure*porosity) - m2);
            jx = jx + Fx;
            m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0) - m4)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
            jy = jy + Fy;
            m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0) - m6)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
            jz = jz + Fz;
            m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0) - m8)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
            m9 = m9 + rlx_setA*((rho0*(2*ux*ux-uy*uy-uz*uz)/porosity) - m9);
            m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10);
            m11 = m11 + rlx_setA*((rho0*(uy*uy-uz*uz)/porosity) - m11);
            m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12);
            m13 = m13 + rlx_setA*((rho0*ux*uy/porosity) - m13);
            m14 = m14 + rlx_setA*((rho0*uy*uz/porosity) - m14);
            m15 = m15 + rlx_setA*((rho0*ux*uz/porosity) - m15);
            m16 = m16 + rlx_setB*( - m16);
            m17 = m17 + rlx_setB*( - m17);
            m18 = m18 + rlx_setB*( - m18);
            //.......................................................................................................

            //.................inverse transformation......................................................
            // q=0
            fq = mrt_V1*rho0-mrt_V2*m1+mrt_V3*m2;
            dist[n] = fq;

            // q = 1
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10);
            dist[1*Np+n] = fq;

            // q=2
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10);
            dist[2*Np+n] = fq;

            // q = 3
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
            dist[3*Np+n] = fq;

            // q = 4
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
            dist[4*Np+n] = fq;

            // q = 5
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
            dist[5*Np+n] = fq;

            // q = 6
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
            dist[6*Np+n] = fq;

            // q = 7
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m16-m17);
            dist[7*Np+n] = fq;

            // q = 8
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m17-m16);
            dist[8*Np+n] = fq;

            // q = 9
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13+0.125*(m16+m17);
            dist[9*Np+n] = fq;

            // q = 10
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13-0.125*(m16+m17);
            dist[10*Np+n] = fq;

            // q = 11
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m18-m16);
            dist[11*Np+n] = fq;

            // q = 12
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m16-m18);
            dist[12*Np+n] = fq;

            // q = 13
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12-0.25*m15-0.125*(m16+m18);
            dist[13*Np+n] = fq;

            // q= 14
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12-0.25*m15+0.125*(m16+m18);
            dist[14*Np+n] = fq;

            // q = 15
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18);
            dist[15*Np+n] = fq;

            // q = 16
            fq =  mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17);
            dist[16*Np+n] = fq;

            // q = 17
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)-mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18);
            dist[17*Np+n] = fq;

            // q = 18
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)-mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18);
            dist[18*Np+n] = fq;
            //........................................................................

            //Update velocity on device
            Velocity[0*Np+n] = ux;
            Velocity[1*Np+n] = uy;
            Velocity[2*Np+n] = uz;
            //Update pressure on device
            Pressure[n] = pressure;

            //-----------------------Mass transport------------------------//
            // Calculate the color gradient
            nx = (2*nB*nA_gradx-2*nA*nB_gradx)/(nA+nB)/(nA+nB); 
            ny = (2*nB*nA_grady-2*nA*nB_grady)/(nA+nB)/(nA+nB); 
            nz = (2*nB*nA_gradz-2*nA*nB_gradz)/(nA+nB)/(nA+nB); 
			//...........Normalize the Color Gradient.................................
			C = sqrt(nx*nx+ny*ny+nz*nz);
			double ColorMag = C;
			if (C==0.0) ColorMag=1.0;
			nx = nx/ColorMag;
			ny = ny/ColorMag;
			nz = nz/ColorMag;		
			if (C == 0.0)	nx = ny = nz = 0.0;

			// Instantiate mass transport distributions
			// Stationary value - distribution 0
			nAB = 1.0/(nA+nB);
			Aq[n] = 0.3333333333333333*nA;
			Bq[n] = 0.3333333333333333*nB;

			//...............................................
			// q = 0,2,4
			// Cq = {1,0,0}, {0,1,0}, {0,0,1}
			delta = beta*nA*nB*nAB*0.1111111111111111*nx;
			if (!(nA*nB*nAB>0)) delta=0;
			a1 = nA*(0.1111111111111111*(1+4.5*ux))+delta;
			b1 = nB*(0.1111111111111111*(1+4.5*ux))-delta;
			a2 = nA*(0.1111111111111111*(1-4.5*ux))-delta;
			b2 = nB*(0.1111111111111111*(1-4.5*ux))+delta;

			Aq[1*Np+n] = a1;
			Bq[1*Np+n] = b1;
			Aq[2*Np+n] = a2;
			Bq[2*Np+n] = b2;

			//...............................................
			// q = 2
			// Cq = {0,1,0}
			delta = beta*nA*nB*nAB*0.1111111111111111*ny;
			if (!(nA*nB*nAB>0)) delta=0;
			a1 = nA*(0.1111111111111111*(1+4.5*uy))+delta;
			b1 = nB*(0.1111111111111111*(1+4.5*uy))-delta;
			a2 = nA*(0.1111111111111111*(1-4.5*uy))-delta;
			b2 = nB*(0.1111111111111111*(1-4.5*uy))+delta;

			Aq[3*Np+n] = a1;
			Bq[3*Np+n] = b1;
			Aq[4*Np+n] = a2;
			Bq[4*Np+n] = b2;
			//...............................................
			// q = 4
			// Cq = {0,0,1}
			delta = beta*nA*nB*nAB*0.1111111111111111*nz;
			if (!(nA*nB*nAB>0)) delta=0;
			a1 = nA*(0.1111111111111111*(1+4.5*uz))+delta;
			b1 = nB*(0.1111111111111111*(1+4.5*uz))-delta;
			a2 = nA*(0.1111111111111111*(1-4.5*uz))-delta;
			b2 = nB*(0.1111111111111111*(1-4.5*uz))+delta;

			Aq[5*Np+n] = a1;
			Bq[5*Np+n] = b1;
			Aq[6*Np+n] = a2;
			Bq[6*Np+n] = b2;
			//...............................................

		}
	}
}

__global__ void dvc_ScaLBL_D3Q19_AAodd_GreyscaleColor(int *neighborList, double *dist, double *Aq, double *Bq, double *Den,
                double *DenGradA, double *DenGradB, double *SolidForce, int start, int finish, int Np,
                double tauA,double tauB,double tauA_eff,double tauB_eff,double rhoA,double rhoB,double Gsc, double Gx, double Gy, double Gz,
                double *Poros,double *Perm, double *Velocity,double *Pressure){

	int n, nread, nr1,nr2,nr3,nr4,nr5,nr6;
	double vx,vy,vz,v_mag;
    double ux,uy,uz,u_mag;
    double pressure;//defined for this incompressible model
	// conserved momemnts
	double jx,jy,jz;
	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
    double fq;
    // currently disable 'GeoFun'
    double GeoFun=0.0;//geometric function from Guo's PRE 66, 036304 (2002)
    double porosity;
    double perm;//voxel permeability
    double c0, c1; //Guo's model parameters
    double Fx, Fy, Fz;//The total body force including Brinkman force and user-specified (Gx,Gy,Gz)
	double tau,tau_eff,rlx_setA,rlx_setB;
    double mu_eff;//effective kinematic viscosity for Darcy term
    double rho0;
    double phi;
    double nx,ny,nz,C;
    double nA,nB;
	double a1,b1,a2,b2,nAB,delta;
    double beta=0.95;
    double nA_gradx,nA_grady,nA_gradz;
    double nB_gradx,nB_grady,nB_gradz;
    double Gff_x,Gff_y,Gff_z;
    double Gfs_x,Gfs_y,Gfs_z;

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

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
	    //........Get 1-D index for this thread....................
	    n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

		if ( n<finish ){		

			// read the component number densities
			nA = Den[n];
			nB = Den[Np + n];
			nA_gradx = DenGradA[n+0*Np];
			nA_grady = DenGradA[n+1*Np];
			nA_gradz = DenGradA[n+2*Np];
			nB_gradx = DenGradB[n+0*Np];
			nB_grady = DenGradB[n+1*Np];
			nB_gradz = DenGradB[n+2*Np];

			// compute phase indicator field
			phi=(nA-nB)/(nA+nB);

			// local density
			rho0=rhoA + 0.5*(1.0-phi)*(rhoB-rhoA);
			// local relaxation time
			tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);
			rlx_setA = 1.f/tau;
			rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
			tau_eff=tauA_eff + 0.5*(1.0-phi)*(tauB_eff-tauA_eff);
            mu_eff = (tau_eff-0.5)/3.f;//kinematic viscosity

            //........................................................................
            //					READ THE DISTRIBUTIONS
            //		(read from opposite array due to previous swap operation)
            //........................................................................
            // q=0
            fq = dist[n];
            m1  = -30.0*fq;
            m2  = 12.0*fq;

            // q=1
            nr1 = neighborList[n]; // neighbor 2 ( > 10Np => odd part of dist)
            fq = dist[nr1]; // reading the f1 data into register fq
            pressure = fq;
            m1 -= 11.0*fq;
            m2 -= 4.0*fq;
            jx = fq;
            m4 = -4.0*fq;
            m9 = 2.0*fq;
            m10 = -4.0*fq;

            // q=2
            nr2 = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
            fq = dist[nr2];  // reading the f2 data into register fq
            pressure += fq;
            m1 -= 11.0*(fq);
            m2 -= 4.0*(fq);
            jx -= fq;
            m4 += 4.0*(fq);
            m9 += 2.0*(fq);
            m10 -= 4.0*(fq);

            // q=3
            nr3 = neighborList[n+2*Np]; // neighbor 4
            fq = dist[nr3];
            pressure += fq;
            m1 -= 11.0*fq;
            m2 -= 4.0*fq;
            jy = fq;
            m6 = -4.0*fq;
            m9 -= fq;
            m10 += 2.0*fq;
            m11 = fq;
            m12 = -2.0*fq;

            // q = 4
            nr4 = neighborList[n+3*Np]; // neighbor 3
            fq = dist[nr4];
            pressure += fq;
            m1 -= 11.0*fq;
            m2 -= 4.0*fq;
            jy -= fq;
            m6 += 4.0*fq;
            m9 -= fq;
            m10 += 2.0*fq;
            m11 += fq;
            m12 -= 2.0*fq;

            // q=5
            nr5 = neighborList[n+4*Np];
            fq = dist[nr5];
            pressure += fq;
            m1 -= 11.0*fq;
            m2 -= 4.0*fq;
            jz = fq;
            m8 = -4.0*fq;
            m9 -= fq;
            m10 += 2.0*fq;
            m11 -= fq;
            m12 += 2.0*fq;

            // q = 6
            nr6 = neighborList[n+5*Np];
            fq = dist[nr6];
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            nread = neighborList[n+16*Np];
            fq = dist[nread];
            pressure += fq;
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
            pressure += fq;
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
            //---------------------------------------------------------------------//

            //---------------- Calculate SC fluid-fluid and fluid-solid forces ---------------//
            // fluid-fluid force
//            Gff_x = -Gsc*nA*nB_gradx*int(phi>0.0)-Gsc*nB*nA_gradx*int(phi<0.0);
//            Gff_y = -Gsc*nA*nB_grady*int(phi>0.0)-Gsc*nB*nA_grady*int(phi<0.0);
//            Gff_z = -Gsc*nA*nB_gradz*int(phi>0.0)-Gsc*nB*nA_gradz*int(phi<0.0);
            Gff_x = -Gsc*(nA*nB_gradx+nB*nA_gradx);
            Gff_y = -Gsc*(nA*nB_grady+nB*nA_grady);
            Gff_z = -Gsc*(nA*nB_gradz+nB*nA_gradz);
            // fluid-solid force
            Gfs_x = (nA-nB)*SolidForce[n+0*Np];    
            Gfs_y = (nA-nB)*SolidForce[n+1*Np];    
            Gfs_z = (nA-nB)*SolidForce[n+2*Np];    

            porosity = Poros[n];
            // use local saturation as an estimation of effective relperm values
            perm = Perm[n]*nA/(nA+nB)*int(phi>0.0)+Perm[n]*nB/(nA+nB)*int(phi<0.0);

            c0 = 0.5*(1.0+porosity*0.5*mu_eff/perm);
            if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
            //GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
            c1 = porosity*0.5*GeoFun/sqrt(perm);
            if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes

            vx = jx/rho0+0.5*(porosity*Gx+Gff_x+Gfs_x);
            vy = jy/rho0+0.5*(porosity*Gy+Gff_y+Gfs_y);
            vz = jz/rho0+0.5*(porosity*Gz+Gff_z+Gfs_z);
            v_mag=sqrt(vx*vx+vy*vy+vz*vz);
            ux = vx/(c0+sqrt(c0*c0+c1*v_mag));
            uy = vy/(c0+sqrt(c0*c0+c1*v_mag));
            uz = vz/(c0+sqrt(c0*c0+c1*v_mag));
            u_mag=sqrt(ux*ux+uy*uy+uz*uz);

            //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
            Fx = rho0*(-porosity*mu_eff/perm*ux - porosity*GeoFun/sqrt(perm)*u_mag*ux + porosity*Gx + Gff_x + Gfs_x);
            Fy = rho0*(-porosity*mu_eff/perm*uy - porosity*GeoFun/sqrt(perm)*u_mag*uy + porosity*Gy + Gff_y + Gfs_y);
            Fz = rho0*(-porosity*mu_eff/perm*uz - porosity*GeoFun/sqrt(perm)*u_mag*uz + porosity*Gz + Gff_z + Gfs_z);
            if (porosity==1.0){
                Fx=rho0*(Gx + Gff_x + Gfs_x);
                Fy=rho0*(Gy + Gff_y + Gfs_y);
                Fz=rho0*(Gz + Gff_z + Gfs_z);
            }

            //Calculate pressure for Incompressible-MRT model
            pressure=0.5/porosity*(pressure-0.5*rho0*u_mag*u_mag/porosity);

//            //..............carry out relaxation process...............................................
//            m1 = m1 + rlx_setA*((-30*rho0+19*(ux*ux+uy*uy+uz*uz)/porosity + 57*pressure*porosity) - m1) 
//                    + (1-0.5*rlx_setA)*38*(Fx*ux+Fy*uy+Fz*uz)/porosity;
//            m2 = m2 + rlx_setA*((12*rho0 - 5.5*(ux*ux+uy*uy+uz*uz)/porosity-27*pressure*porosity) - m2)
//                    + (1-0.5*rlx_setA)*11*(-Fx*ux-Fy*uy-Fz*uz)/porosity;
//            jx = jx + Fx;
//            m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0) - m4)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
//            jy = jy + Fy;
//            m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0) - m6)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
//            jz = jz + Fz;
//            m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0) - m8)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
//            m9 = m9 + rlx_setA*((rho0*(2*ux*ux-uy*uy-uz*uz)/porosity) - m9)
//                    + (1-0.5*rlx_setA)*(4*Fx*ux-2*Fy*uy-2*Fz*uz)/porosity;
//            m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10)
//                      + (1-0.5*rlx_setA)*(-2*Fx*ux+Fy*uy+Fz*uz)/porosity;
//            m11 = m11 + rlx_setA*((rho0*(uy*uy-uz*uz)/porosity) - m11)
//                      + (1-0.5*rlx_setA)*(2*Fy*uy-2*Fz*uz)/porosity;
//            m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12)
//                      + (1-0.5*rlx_setA)*(-Fy*uy+Fz*uz)/porosity;
//            m13 = m13 + rlx_setA*((rho0*ux*uy/porosity) - m13)
//                      + (1-0.5*rlx_setA)*(Fy*ux+Fx*uy)/porosity;
//            m14 = m14 + rlx_setA*((rho0*uy*uz/porosity) - m14)
//                      + (1-0.5*rlx_setA)*(Fz*uy+Fy*uz)/porosity;
//            m15 = m15 + rlx_setA*((rho0*ux*uz/porosity) - m15)
//                      + (1-0.5*rlx_setA)*(Fz*ux+Fx*uz)/porosity;
//            m16 = m16 + rlx_setB*( - m16);
//            m17 = m17 + rlx_setB*( - m17);
//            m18 = m18 + rlx_setB*( - m18);
//            //.......................................................................................................
           
            //-------------------- IMRT collison where body force has NO higher-order terms -------------//
            //..............carry out relaxation process...............................................
            m1 = m1 + rlx_setA*((-30*rho0+19*(ux*ux+uy*uy+uz*uz)/porosity + 57*pressure*porosity) - m1);
            m2 = m2 + rlx_setA*((12*rho0 - 5.5*(ux*ux+uy*uy+uz*uz)/porosity-27*pressure*porosity) - m2);
            jx = jx + Fx;
            m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0) - m4)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
            jy = jy + Fy;
            m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0) - m6)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
            jz = jz + Fz;
            m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0) - m8)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
            m9 = m9 + rlx_setA*((rho0*(2*ux*ux-uy*uy-uz*uz)/porosity) - m9);
            m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10);
            m11 = m11 + rlx_setA*((rho0*(uy*uy-uz*uz)/porosity) - m11);
            m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12);
            m13 = m13 + rlx_setA*((rho0*ux*uy/porosity) - m13);
            m14 = m14 + rlx_setA*((rho0*uy*uz/porosity) - m14);
            m15 = m15 + rlx_setA*((rho0*ux*uz/porosity) - m15);
            m16 = m16 + rlx_setB*( - m16);
            m17 = m17 + rlx_setB*( - m17);
            m18 = m18 + rlx_setB*( - m18);
            //.......................................................................................................


            //.................inverse transformation......................................................
            // q=0
            fq = mrt_V1*rho0-mrt_V2*m1+mrt_V3*m2;
            dist[n] = fq;

            // q = 1
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10);
            //nread = neighborList[n+Np];
            dist[nr2] = fq;

            // q=2
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10);
            //nread = neighborList[n];
            dist[nr1] = fq;

            // q = 3
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
            //nread = neighborList[n+3*Np];
            dist[nr4] = fq;

            // q = 4
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
            //nread = neighborList[n+2*Np];
            dist[nr3] = fq;

            // q = 5
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
            //nread = neighborList[n+5*Np];
            dist[nr6] = fq;

            // q = 6
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
            //nread = neighborList[n+4*Np];
            dist[nr5] = fq;

            // q = 7
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m16-m17);
            nread = neighborList[n+7*Np];
            dist[nread] = fq;

            // q = 8
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m17-m16);
            nread = neighborList[n+6*Np];
            dist[nread] = fq;

            // q = 9
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13+0.125*(m16+m17);
            nread = neighborList[n+9*Np];
            dist[nread] = fq;

            // q = 10
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13-0.125*(m16+m17);
            nread = neighborList[n+8*Np];
            dist[nread] = fq;

            // q = 11
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m18-m16);
            nread = neighborList[n+11*Np];
            dist[nread] = fq;

            // q = 12
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m16-m18);
            nread = neighborList[n+10*Np];
            dist[nread]= fq;

            // q = 13
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12-0.25*m15-0.125*(m16+m18);
            nread = neighborList[n+13*Np];
            dist[nread] = fq;

            // q= 14
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12-0.25*m15+0.125*(m16+m18);
            nread = neighborList[n+12*Np];
            dist[nread] = fq;

            // q = 15
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18);
            nread = neighborList[n+15*Np];
            dist[nread] = fq;

            // q = 16
            fq =  mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17);
            nread = neighborList[n+14*Np];
            dist[nread] = fq;

            // q = 17
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)-mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18);
            nread = neighborList[n+17*Np];
            dist[nread] = fq;

            // q = 18
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)-mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18);
            nread = neighborList[n+16*Np];
            dist[nread] = fq;
            //........................................................................

            //Update velocity on device
            Velocity[0*Np+n] = ux;
            Velocity[1*Np+n] = uy;
            Velocity[2*Np+n] = uz;
            //Update pressure on device
            Pressure[n] = pressure;

            //-----------------------Mass transport------------------------//
            // Calculate the color gradient
            nx = (2*nB*nA_gradx-2*nA*nB_gradx)/(nA+nB)/(nA+nB); 
            ny = (2*nB*nA_grady-2*nA*nB_grady)/(nA+nB)/(nA+nB); 
            nz = (2*nB*nA_gradz-2*nA*nB_gradz)/(nA+nB)/(nA+nB); 
			//...........Normalize the Color Gradient.................................
			C = sqrt(nx*nx+ny*ny+nz*nz);
			double ColorMag = C;
			if (C==0.0) ColorMag=1.0;
			nx = nx/ColorMag;
			ny = ny/ColorMag;
			nz = nz/ColorMag;		
			if (C == 0.0)	nx = ny = nz = 0.0;

			// Instantiate mass transport distributions
			// Stationary value - distribution 0
			nAB = 1.0/(nA+nB);
			Aq[n] = 0.3333333333333333*nA;
			Bq[n] = 0.3333333333333333*nB;

			//...............................................
			// q = 0,2,4
			// Cq = {1,0,0}, {0,1,0}, {0,0,1}
			delta = beta*nA*nB*nAB*0.1111111111111111*nx;
			if (!(nA*nB*nAB>0)) delta=0;
			a1 = nA*(0.1111111111111111*(1+4.5*ux))+delta;
			b1 = nB*(0.1111111111111111*(1+4.5*ux))-delta;
			a2 = nA*(0.1111111111111111*(1-4.5*ux))-delta;
			b2 = nB*(0.1111111111111111*(1-4.5*ux))+delta;

			// q = 1
			//nread = neighborList[n+Np];
			Aq[nr2] = a1;
			Bq[nr2] = b1;
			// q=2
			//nread = neighborList[n];
			Aq[nr1] = a2;
			Bq[nr1] = b2;

			//...............................................
			// Cq = {0,1,0}
			delta = beta*nA*nB*nAB*0.1111111111111111*ny;
			if (!(nA*nB*nAB>0)) delta=0;
			a1 = nA*(0.1111111111111111*(1+4.5*uy))+delta;
			b1 = nB*(0.1111111111111111*(1+4.5*uy))-delta;
			a2 = nA*(0.1111111111111111*(1-4.5*uy))-delta;
			b2 = nB*(0.1111111111111111*(1-4.5*uy))+delta;

			// q = 3
			//nread = neighborList[n+3*Np];
			Aq[nr4] = a1;
			Bq[nr4] = b1;
			// q = 4
			//nread = neighborList[n+2*Np];
			Aq[nr3] = a2;
			Bq[nr3] = b2;

			//...............................................
			// q = 4
			// Cq = {0,0,1}
			delta = beta*nA*nB*nAB*0.1111111111111111*nz;
			if (!(nA*nB*nAB>0)) delta=0;
			a1 = nA*(0.1111111111111111*(1+4.5*uz))+delta;
			b1 = nB*(0.1111111111111111*(1+4.5*uz))-delta;
			a2 = nA*(0.1111111111111111*(1-4.5*uz))-delta;
			b2 = nB*(0.1111111111111111*(1-4.5*uz))+delta;

			// q = 5
			//nread = neighborList[n+5*Np];
			Aq[nr6] = a1;
			Bq[nr6] = b1;
			// q = 6
			//nread = neighborList[n+4*Np];
			Aq[nr5] = a2;
			Bq[nr5] = b2;
			//...............................................
		}
	}
}

//__global__ void dvc_ScaLBL_D3Q19_AAeven_GreyscaleColorChem(double *dist, double *Aq, double *Bq, double *Den,double *SolidForce, int start, int finish, int Np,
//                double tauA,double tauB,double tauA_eff,double tauB_eff,double rhoA,double rhoB,double gamma,double kappaA,double kappaB,double lambdaA,double lambdaB,
//                double Gx, double Gy, double Gz,
//                double *Poros,double *Perm, double *Velocity,double *Pressure,double *PressureGrad,double *PressTensorGrad,double *PhiLap){
//	int n;
//	double vx,vy,vz,v_mag;
//    double ux,uy,uz,u_mag;
//    double pressure;//defined for this incompressible model
//	// conserved momemnts
//	double jx,jy,jz;
//	// non-conserved moments
//	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
//    double fq;
//    // currently disable 'GeoFun'
//    double GeoFun=0.0;//geometric function from Guo's PRE 66, 036304 (2002)
//    double porosity;
//    double perm;//voxel permeability
//    double c0, c1; //Guo's model parameters
//    double Fx, Fy, Fz;//The total body force including Brinkman force and user-specified (Gx,Gy,Gz)
//	double tau,tau_eff,rlx_setA,rlx_setB;
//    double mu_eff;//effective kinematic viscosity for Darcy term
//    double rho0;
//    double phi;
//    double phi_lap;//laplacian of phase field
//    double nA,nB;
//	double a1,b1,a2,b2;
//    double Gfs_x,Gfs_y,Gfs_z;
//    double Gff_x,Gff_y,Gff_z;
//    double chem_a,chem_b;
//    double rlx_massA,rlx_massB;
//    // *---------------------------------Pressure Tensor Gradient------------------------------------*//
//    double Pxx_x,Pyy_y,Pzz_z;
//    double Pxy_x,Pxy_y;
//    double Pyz_y,Pyz_z;
//    double Pxz_x,Pxz_z;
//    double px,py,pz; //pressure gradient
//
//
//	const double mrt_V1=0.05263157894736842;
//	const double mrt_V2=0.012531328320802;
//	const double mrt_V3=0.04761904761904762;
//	const double mrt_V4=0.004594820384294068;
//	const double mrt_V5=0.01587301587301587;
//	const double mrt_V6=0.0555555555555555555555555;
//	const double mrt_V7=0.02777777777777778;
//	const double mrt_V8=0.08333333333333333;
//	const double mrt_V9=0.003341687552213868;
//	const double mrt_V10=0.003968253968253968;
//	const double mrt_V11=0.01388888888888889;
//	const double mrt_V12=0.04166666666666666;
//
//
//	int S = Np/NBLOCKS/NTHREADS + 1;
//	for (int s=0; s<S; s++){
//	    //........Get 1-D index for this thread....................
//	    n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
//
//		if ( n<finish ){
//
//			// read the component number densities
//			nA = Den[n];
//			nB = Den[Np + n];
//			// compute phase indicator field
//			phi=(nA-nB)/(nA+nB);
//            // load laplacian of phase field
//            phi_lap = PhiLap[n];
//            // Load voxel porosity and perm
//            porosity = Poros[n];
//            // use local saturation as an estimation of effective relperm values
//            perm = Perm[n]*nA/(nA+nB)*int(phi>0.0)+Perm[n]*nB/(nA+nB)*int(phi<0.0);
//
//            //Load pressure gradient
//            px=PressureGrad[0*Np+n];
//            py=PressureGrad[1*Np+n];
//            pz=PressureGrad[2*Np+n];
//
//            //Load pressure tensor gradient
//            //For reference full list of PressTensorGrad
//            //PressTensorGrad[n+0*Np]  = Pxx_x
//            //PressTensorGrad[n+1*Np]  = Pxx_y
//            //PressTensorGrad[n+2*Np]  = Pxx_z
//            //PressTensorGrad[n+3*Np]  = Pyy_x
//            //PressTensorGrad[n+4*Np]  = Pyy_y
//            //PressTensorGrad[n+5*Np]  = Pyy_z
//            //PressTensorGrad[n+6*Np]  = Pzz_x
//            //PressTensorGrad[n+7*Np]  = Pzz_y
//            //PressTensorGrad[n+8*Np]  = Pzz_z
//            //PressTensorGrad[n+9*Np]  = Pxy_x
//            //PressTensorGrad[n+10*Np] = Pxy_y
//            //PressTensorGrad[n+11*Np] = Pxy_z
//            //PressTensorGrad[n+12*Np] = Pyz_x
//            //PressTensorGrad[n+13*Np] = Pyz_y
//            //PressTensorGrad[n+14*Np] = Pyz_z
//            //PressTensorGrad[n+15*Np] = Pxz_x
//            //PressTensorGrad[n+16*Np] = Pxz_y
//            //PressTensorGrad[n+17*Np] = Pxz_z
//            Pxx_x = PressTensorGrad[0*Np+n];
//            Pyy_y = PressTensorGrad[4*Np+n];
//            Pzz_z = PressTensorGrad[8*Np+n];
//            Pxy_x = PressTensorGrad[9*Np+n];
//            Pxz_x = PressTensorGrad[15*Np+n];
//		    Pxy_y = PressTensorGrad[10*Np+n];
//		    Pyz_y = PressTensorGrad[13*Np+n];
//		    Pyz_z = PressTensorGrad[14*Np+n];
//		    Pxz_z = PressTensorGrad[17*Np+n];
//		    //............Compute the fluid-fluid force (gfx,gfy,gfz)...................................
//            //TODO double check if you need porosity as a fre-factor
//            Gff_x = porosity*px-(Pxx_x+Pxy_y+Pxz_z);
//            Gff_y = porosity*py-(Pxy_x+Pyy_y+Pyz_z);
//            Gff_z = porosity*pz-(Pxz_x+Pyz_y+Pzz_z);
//            // fluid-solid force
//            Gfs_x = (nA-nB)*SolidForce[n+0*Np];    
//            Gfs_y = (nA-nB)*SolidForce[n+1*Np];    
//            Gfs_z = (nA-nB)*SolidForce[n+2*Np];    
//
//			// local density
//			rho0=rhoA + 0.5*(1.0-phi)*(rhoB-rhoA);
//			// local relaxation time
//			tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);
//			rlx_setA = 1.f/tau;
//			rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
//			tau_eff=tauA_eff + 0.5*(1.0-phi)*(tauB_eff-tauA_eff);
//            mu_eff = (tau_eff-0.5)/3.f;//kinematic viscosity
//
//
//            //........................................................................
//            //					READ THE DISTRIBUTIONS
//            //		(read from opposite array due to previous swap operation)
//            //........................................................................
//            // q=0
//            fq = dist[n];
//            m1  = -30.0*fq;
//            m2  = 12.0*fq;
//
//            // q=1
//            fq = dist[2*Np+n];
//            pressure = fq;
//            m1 -= 11.0*fq;
//            m2 -= 4.0*fq;
//            jx = fq;
//            m4 = -4.0*fq;
//            m9 = 2.0*fq;
//            m10 = -4.0*fq;
//
//            // f2 = dist[10*Np+n];
//            fq = dist[1*Np+n];
//            pressure += fq;
//            m1 -= 11.0*(fq);
//            m2 -= 4.0*(fq);
//            jx -= fq;
//            m4 += 4.0*(fq);
//            m9 += 2.0*(fq);
//            m10 -= 4.0*(fq);
//
//            // q=3
//            fq = dist[4*Np+n];
//            pressure += fq;
//            m1 -= 11.0*fq;
//            m2 -= 4.0*fq;
//            jy = fq;
//            m6 = -4.0*fq;
//            m9 -= fq;
//            m10 += 2.0*fq;
//            m11 = fq;
//            m12 = -2.0*fq;
//
//            // q = 4
//            fq = dist[3*Np+n];
//            pressure += fq;
//            m1 -= 11.0*fq;
//            m2 -= 4.0*fq;
//            jy -= fq;
//            m6 += 4.0*fq;
//            m9 -= fq;
//            m10 += 2.0*fq;
//            m11 += fq;
//            m12 -= 2.0*fq;
//
//            // q=5
//            fq = dist[6*Np+n];
//            pressure += fq;
//            m1 -= 11.0*fq;
//            m2 -= 4.0*fq;
//            jz = fq;
//            m8 = -4.0*fq;
//            m9 -= fq;
//            m10 += 2.0*fq;
//            m11 -= fq;
//            m12 += 2.0*fq;
//
//            // q = 6
//            fq = dist[5*Np+n];
//            pressure += fq;
//            m1 -= 11.0*fq;
//            m2 -= 4.0*fq;
//            jz -= fq;
//            m8 += 4.0*fq;
//            m9 -= fq;
//            m10 += 2.0*fq;
//            m11 -= fq;
//            m12 += 2.0*fq;
//
//            // q=7
//            fq = dist[8*Np+n];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jx += fq;
//            m4 += fq;
//            jy += fq;
//            m6 += fq;
//            m9  += fq;
//            m10 += fq;
//            m11 += fq;
//            m12 += fq;
//            m13 = fq;
//            m16 = fq;
//            m17 = -fq;
//
//            // q = 8
//            fq = dist[7*Np+n];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jx -= fq;
//            m4 -= fq;
//            jy -= fq;
//            m6 -= fq;
//            m9 += fq;
//            m10 += fq;
//            m11 += fq;
//            m12 += fq;
//            m13 += fq;
//            m16 -= fq;
//            m17 += fq;
//
//            // q=9
//            fq = dist[10*Np+n];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jx += fq;
//            m4 += fq;
//            jy -= fq;
//            m6 -= fq;
//            m9 += fq;
//            m10 += fq;
//            m11 += fq;
//            m12 += fq;
//            m13 -= fq;
//            m16 += fq;
//            m17 += fq;
//
//            // q = 10
//            fq = dist[9*Np+n];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jx -= fq;
//            m4 -= fq;
//            jy += fq;
//            m6 += fq;
//            m9 += fq;
//            m10 += fq;
//            m11 += fq;
//            m12 += fq;
//            m13 -= fq;
//            m16 -= fq;
//            m17 -= fq;
//
//            // q=11
//            fq = dist[12*Np+n];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jx += fq;
//            m4 += fq;
//            jz += fq;
//            m8 += fq;
//            m9 += fq;
//            m10 += fq;
//            m11 -= fq;
//            m12 -= fq;
//            m15 = fq;
//            m16 -= fq;
//            m18 = fq;
//
//            // q=12
//            fq = dist[11*Np+n];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jx -= fq;
//            m4 -= fq;
//            jz -= fq;
//            m8 -= fq;
//            m9 += fq;
//            m10 += fq;
//            m11 -= fq;
//            m12 -= fq;
//            m15 += fq;
//            m16 += fq;
//            m18 -= fq;
//
//            // q=13
//            fq = dist[14*Np+n];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jx += fq;
//            m4 += fq;
//            jz -= fq;
//            m8 -= fq;
//            m9 += fq;
//            m10 += fq;
//            m11 -= fq;
//            m12 -= fq;
//            m15 -= fq;
//            m16 -= fq;
//            m18 -= fq;
//
//            // q=14
//            fq = dist[13*Np+n];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jx -= fq;
//            m4 -= fq;
//            jz += fq;
//            m8 += fq;
//            m9 += fq;
//            m10 += fq;
//            m11 -= fq;
//            m12 -= fq;
//            m15 -= fq;
//            m16 += fq;
//            m18 += fq;
//
//            // q=15
//            fq = dist[16*Np+n];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jy += fq;
//            m6 += fq;
//            jz += fq;
//            m8 += fq;
//            m9 -= 2.0*fq;
//            m10 -= 2.0*fq;
//            m14 = fq;
//            m17 += fq;
//            m18 -= fq;
//
//            // q=16
//            fq = dist[15*Np+n];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jy -= fq;
//            m6 -= fq;
//            jz -= fq;
//            m8 -= fq;
//            m9 -= 2.0*fq;
//            m10 -= 2.0*fq;
//            m14 += fq;
//            m17 -= fq;
//            m18 += fq;
//
//            // q=17
//            fq = dist[18*Np+n];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jy += fq;
//            m6 += fq;
//            jz -= fq;
//            m8 -= fq;
//            m9 -= 2.0*fq;
//            m10 -= 2.0*fq;
//            m14 -= fq;
//            m17 += fq;
//            m18 += fq;
//
//            // q=18
//            fq = dist[17*Np+n];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jy -= fq;
//            m6 -= fq;
//            jz += fq;
//            m8 += fq;
//            m9 -= 2.0*fq;
//            m10 -= 2.0*fq;
//            m14 -= fq;
//            m17 -= fq;
//            m18 -= fq;
//            //---------------------------------------------------------------------//
//
//            c0 = 0.5*(1.0+porosity*0.5*mu_eff/perm);
//            if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
//            //GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
//            c1 = porosity*0.5*GeoFun/sqrt(perm);
//            if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes
//
//            vx = jx/rho0+0.5*(porosity*Gx+Gff_x+Gfs_x);
//            vy = jy/rho0+0.5*(porosity*Gy+Gff_y+Gfs_y);
//            vz = jz/rho0+0.5*(porosity*Gz+Gff_z+Gfs_z);
//            v_mag=sqrt(vx*vx+vy*vy+vz*vz);
//            ux = vx/(c0+sqrt(c0*c0+c1*v_mag));
//            uy = vy/(c0+sqrt(c0*c0+c1*v_mag));
//            uz = vz/(c0+sqrt(c0*c0+c1*v_mag));
//            u_mag=sqrt(ux*ux+uy*uy+uz*uz);
//
//            //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
//            Fx = rho0*(-porosity*mu_eff/perm*ux - porosity*GeoFun/sqrt(perm)*u_mag*ux + porosity*Gx + Gff_x + Gfs_x);
//            Fy = rho0*(-porosity*mu_eff/perm*uy - porosity*GeoFun/sqrt(perm)*u_mag*uy + porosity*Gy + Gff_y + Gfs_y);
//            Fz = rho0*(-porosity*mu_eff/perm*uz - porosity*GeoFun/sqrt(perm)*u_mag*uz + porosity*Gz + Gff_z + Gfs_z);
//            if (porosity==1.0){
//                Fx=rho0*(Gx + Gff_x + Gfs_x);
//                Fy=rho0*(Gy + Gff_y + Gfs_y);
//                Fz=rho0*(Gz + Gff_z + Gfs_z);
//            }
//
//            //Calculate pressure for Incompressible-MRT model
//            pressure=0.5/porosity*(pressure-0.5*rho0*u_mag*u_mag/porosity);
//
////            //..............carry out relaxation process...............................................
////            m1 = m1 + rlx_setA*((-30*rho0+19*(ux*ux+uy*uy+uz*uz)/porosity + 57*pressure*porosity) - m1) 
////                    + (1-0.5*rlx_setA)*38*(Fx*ux+Fy*uy+Fz*uz)/porosity;
////            m2 = m2 + rlx_setA*((12*rho0 - 5.5*(ux*ux+uy*uy+uz*uz)/porosity-27*pressure*porosity) - m2)
////                    + (1-0.5*rlx_setA)*11*(-Fx*ux-Fy*uy-Fz*uz)/porosity;
////            jx = jx + Fx;
////            m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0) - m4)
////                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
////            jy = jy + Fy;
////            m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0) - m6)
////                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
////            jz = jz + Fz;
////            m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0) - m8)
////                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
////            m9 = m9 + rlx_setA*((rho0*(2*ux*ux-uy*uy-uz*uz)/porosity) - m9)
////                    + (1-0.5*rlx_setA)*(4*Fx*ux-2*Fy*uy-2*Fz*uz)/porosity;
////            m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10)
////                      + (1-0.5*rlx_setA)*(-2*Fx*ux+Fy*uy+Fz*uz)/porosity;
////            m11 = m11 + rlx_setA*((rho0*(uy*uy-uz*uz)/porosity) - m11)
////                      + (1-0.5*rlx_setA)*(2*Fy*uy-2*Fz*uz)/porosity;
////            m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12)
////                      + (1-0.5*rlx_setA)*(-Fy*uy+Fz*uz)/porosity;
////            m13 = m13 + rlx_setA*((rho0*ux*uy/porosity) - m13)
////                      + (1-0.5*rlx_setA)*(Fy*ux+Fx*uy)/porosity;
////            m14 = m14 + rlx_setA*((rho0*uy*uz/porosity) - m14)
////                      + (1-0.5*rlx_setA)*(Fz*uy+Fy*uz)/porosity;
////            m15 = m15 + rlx_setA*((rho0*ux*uz/porosity) - m15)
////                      + (1-0.5*rlx_setA)*(Fz*ux+Fx*uz)/porosity;
////            m16 = m16 + rlx_setB*( - m16);
////            m17 = m17 + rlx_setB*( - m17);
////            m18 = m18 + rlx_setB*( - m18);
////            //.......................................................................................................
//
//            //-------------------- IMRT collison where body force has NO higher-order terms -------------//
//            //..............carry out relaxation process...............................................
//            m1 = m1 + rlx_setA*((-30*rho0+19*(ux*ux+uy*uy+uz*uz)/porosity + 57*pressure*porosity) - m1);
//            m2 = m2 + rlx_setA*((12*rho0 - 5.5*(ux*ux+uy*uy+uz*uz)/porosity-27*pressure*porosity) - m2);
//            jx = jx + Fx;
//            m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0) - m4)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
//            jy = jy + Fy;
//            m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0) - m6)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
//            jz = jz + Fz;
//            m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0) - m8)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
//            m9 = m9 + rlx_setA*((rho0*(2*ux*ux-uy*uy-uz*uz)/porosity) - m9);
//            m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10);
//            m11 = m11 + rlx_setA*((rho0*(uy*uy-uz*uz)/porosity) - m11);
//            m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12);
//            m13 = m13 + rlx_setA*((rho0*ux*uy/porosity) - m13);
//            m14 = m14 + rlx_setA*((rho0*uy*uz/porosity) - m14);
//            m15 = m15 + rlx_setA*((rho0*ux*uz/porosity) - m15);
//            m16 = m16 + rlx_setB*( - m16);
//            m17 = m17 + rlx_setB*( - m17);
//            m18 = m18 + rlx_setB*( - m18);
//            //.......................................................................................................
//
//            //.................inverse transformation......................................................
//            // q=0
//            fq = mrt_V1*rho0-mrt_V2*m1+mrt_V3*m2;
//            dist[n] = fq;
//
//            // q = 1
//            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10);
//            dist[1*Np+n] = fq;
//
//            // q=2
//            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10);
//            dist[2*Np+n] = fq;
//
//            // q = 3
//            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
//            dist[3*Np+n] = fq;
//
//            // q = 4
//            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
//            dist[4*Np+n] = fq;
//
//            // q = 5
//            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
//            dist[5*Np+n] = fq;
//
//            // q = 6
//            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
//            dist[6*Np+n] = fq;
//
//            // q = 7
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m16-m17);
//            dist[7*Np+n] = fq;
//
//            // q = 8
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m17-m16);
//            dist[8*Np+n] = fq;
//
//            // q = 9
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13+0.125*(m16+m17);
//            dist[9*Np+n] = fq;
//
//            // q = 10
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13-0.125*(m16+m17);
//            dist[10*Np+n] = fq;
//
//            // q = 11
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m18-m16);
//            dist[11*Np+n] = fq;
//
//            // q = 12
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m16-m18);
//            dist[12*Np+n] = fq;
//
//            // q = 13
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12-0.25*m15-0.125*(m16+m18);
//            dist[13*Np+n] = fq;
//
//            // q= 14
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12-0.25*m15+0.125*(m16+m18);
//            dist[14*Np+n] = fq;
//
//            // q = 15
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18);
//            dist[15*Np+n] = fq;
//
//            // q = 16
//            fq =  mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17);
//            dist[16*Np+n] = fq;
//
//            // q = 17
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)-mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18);
//            dist[17*Np+n] = fq;
//
//            // q = 18
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)-mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18);
//            dist[18*Np+n] = fq;
//            //........................................................................
//
//            //Update velocity on device
//            Velocity[0*Np+n] = ux;
//            Velocity[1*Np+n] = uy;
//            Velocity[2*Np+n] = uz;
//            //Update pressure on device
//            Pressure[n] = pressure;
//
//            //-----------------------Mass transport------------------------//
//            // calcuale chemical potential
//            chem_a = lambdaA*(nA*nA*nA-1.5*nA*nA+0.5*nA)-0.25*kappaA*phi_lap;
//            chem_b = -lambdaB*(nB*nB*nB-1.5*nB*nB+0.5*nB)-0.25*kappaB*phi_lap;
//            rlx_massA = 3.f-sqrt(3.f);
//            rlx_massB = 3.f-sqrt(3.f);
//
//			//...............................................
//			// q = 0,2,4
//			// Cq = {1,0,0}, {0,1,0}, {0,0,1}
//			a1 = Aq[1*Np+n];
//			b1 = Bq[1*Np+n];
//			a2 = Aq[2*Np+n];
//			b2 = Bq[2*Np+n];
//			a1 = (1.0-rlx_massA)*a1+rlx_massA*(0.1111111111111111*4.5*(gamma*chem_a+nA*ux));
//			b1 = (1.0-rlx_massB)*b1+rlx_massB*(0.1111111111111111*4.5*(gamma*chem_b+nB*ux));
//			a2 = (1.0-rlx_massA)*a2+rlx_massA*(0.1111111111111111*4.5*(gamma*chem_a-nA*ux));
//			b2 = (1.0-rlx_massB)*b2+rlx_massB*(0.1111111111111111*4.5*(gamma*chem_b-nB*ux));
//
//			Aq[1*Np+n] = a1;
//			Bq[1*Np+n] = b1;
//			Aq[2*Np+n] = a2;
//			Bq[2*Np+n] = b2;
//
//			//...............................................
//			// q = 2
//			// Cq = {0,1,0}
//			a1 = Aq[3*Np+n];
//			b1 = Bq[3*Np+n];
//			a2 = Aq[4*Np+n];
//			b2 = Bq[4*Np+n];
//			a1 = (1.0-rlx_massA)*a1+rlx_massA*(0.1111111111111111*4.5*(gamma*chem_a+nA*uy));
//			b1 = (1.0-rlx_massB)*b1+rlx_massB*(0.1111111111111111*4.5*(gamma*chem_b+nB*uy));
//			a2 = (1.0-rlx_massA)*a2+rlx_massA*(0.1111111111111111*4.5*(gamma*chem_a-nA*uy));
//			b2 = (1.0-rlx_massB)*b2+rlx_massB*(0.1111111111111111*4.5*(gamma*chem_b-nB*uy));
//
//			Aq[3*Np+n] = a1;
//			Bq[3*Np+n] = b1;
//			Aq[4*Np+n] = a2;
//			Bq[4*Np+n] = b2;
//			//...............................................
//			// q = 4
//			// Cq = {0,0,1}
//			a1 = Aq[5*Np+n];
//			b1 = Bq[5*Np+n];
//			a2 = Aq[6*Np+n];
//			b2 = Bq[6*Np+n];
//			a1 = (1.0-rlx_massA)*a1+rlx_massA*(0.1111111111111111*4.5*(gamma*chem_a+nA*uz));
//			b1 = (1.0-rlx_massB)*b1+rlx_massB*(0.1111111111111111*4.5*(gamma*chem_b+nB*uz));
//			a2 = (1.0-rlx_massA)*a2+rlx_massA*(0.1111111111111111*4.5*(gamma*chem_a-nA*uz));
//			b2 = (1.0-rlx_massB)*b2+rlx_massB*(0.1111111111111111*4.5*(gamma*chem_b-nB*uz));
//
//			Aq[5*Np+n] = a1;
//			Bq[5*Np+n] = b1;
//			Aq[6*Np+n] = a2;
//			Bq[6*Np+n] = b2;
//			//...............................................
//
//			// Instantiate mass transport distributions
//			// Stationary value - distribution 0
//            a1=Aq[n];
//            b1=Bq[n];
//			Aq[n] = (1.0-rlx_massA)*a1+rlx_massA*(nA-3.0*gamma*chem_a);
//			Bq[n] = (1.0-rlx_massB)*b1+rlx_massB*(nB-3.0*gamma*chem_b);
//
//
//		}
//	}
//}

//__global__ void dvc_ScaLBL_D3Q19_AAodd_GreyscaleColorChem(int *neighborList, double *dist, double *Aq, double *Bq, double *Den,double *SolidForce, int start, int finish, int Np,
//                double tauA,double tauB,double tauA_eff,double tauB_eff,double rhoA,double rhoB,double gamma,double kappaA,double kappaB,double lambdaA,double lambdaB,
//                double Gx, double Gy, double Gz,
//                double *Poros,double *Perm, double *Velocity,double *Pressure,double *PressureGrad,double *PressTensorGrad,double *PhiLap){
//
//	int n, nread, nr1,nr2,nr3,nr4,nr5,nr6;
//	double vx,vy,vz,v_mag;
//    double ux,uy,uz,u_mag;
//    double pressure;//defined for this incompressible model
//	// conserved momemnts
//	double jx,jy,jz;
//	// non-conserved moments
//	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
//    double fq;
//    // currently disable 'GeoFun'
//    double GeoFun=0.0;//geometric function from Guo's PRE 66, 036304 (2002)
//    double porosity;
//    double perm;//voxel permeability
//    double c0, c1; //Guo's model parameters
//    double Fx, Fy, Fz;//The total body force including Brinkman force and user-specified (Gx,Gy,Gz)
//	double tau,tau_eff,rlx_setA,rlx_setB;
//    double mu_eff;//effective kinematic viscosity for Darcy term
//    double rho0;
//    double phi;
//    double phi_lap;//laplacian of phase field
//    double nA,nB;
//	double a1,b1,a2,b2;
//    double Gfs_x,Gfs_y,Gfs_z;
//    double Gff_x,Gff_y,Gff_z;
//    double chem_a,chem_b;
//    double rlx_massA,rlx_massB;
//    // *---------------------------------Pressure Tensor Gradient------------------------------------*//
//    double Pxx_x,Pyy_y,Pzz_z;
//    double Pxy_x,Pxy_y;
//    double Pyz_y,Pyz_z;
//    double Pxz_x,Pxz_z;
//    double px,py,pz; //pressure gradient
//
//	const double mrt_V1=0.05263157894736842;
//	const double mrt_V2=0.012531328320802;
//	const double mrt_V3=0.04761904761904762;
//	const double mrt_V4=0.004594820384294068;
//	const double mrt_V5=0.01587301587301587;
//	const double mrt_V6=0.0555555555555555555555555;
//	const double mrt_V7=0.02777777777777778;
//	const double mrt_V8=0.08333333333333333;
//	const double mrt_V9=0.003341687552213868;
//	const double mrt_V10=0.003968253968253968;
//	const double mrt_V11=0.01388888888888889;
//	const double mrt_V12=0.04166666666666666;
//
//	int S = Np/NBLOCKS/NTHREADS + 1;
//	for (int s=0; s<S; s++){
//	    //........Get 1-D index for this thread....................
//	    n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
//
//		if ( n<finish ){		
//
//			// read the component number densities
//			nA = Den[n];
//			nB = Den[Np + n];
//			// compute phase indicator field
//			phi=(nA-nB)/(nA+nB);
//            // load laplacian of phase field
//            phi_lap = PhiLap[n];
//            // Load voxel porosity and perm
//            porosity = Poros[n];
//            // use local saturation as an estimation of effective relperm values
//            perm = Perm[n]*nA/(nA+nB)*int(phi>0.0)+Perm[n]*nB/(nA+nB)*int(phi<0.0);
//
//            //Load pressure gradient
//            px=PressureGrad[0*Np+n];
//            py=PressureGrad[1*Np+n];
//            pz=PressureGrad[2*Np+n];
//
//            //Load pressure tensor gradient
//            //For reference full list of PressTensorGrad
//            //PressTensorGrad[n+0*Np]  = Pxx_x
//            //PressTensorGrad[n+1*Np]  = Pxx_y
//            //PressTensorGrad[n+2*Np]  = Pxx_z
//            //PressTensorGrad[n+3*Np]  = Pyy_x
//            //PressTensorGrad[n+4*Np]  = Pyy_y
//            //PressTensorGrad[n+5*Np]  = Pyy_z
//            //PressTensorGrad[n+6*Np]  = Pzz_x
//            //PressTensorGrad[n+7*Np]  = Pzz_y
//            //PressTensorGrad[n+8*Np]  = Pzz_z
//            //PressTensorGrad[n+9*Np]  = Pxy_x
//            //PressTensorGrad[n+10*Np] = Pxy_y
//            //PressTensorGrad[n+11*Np] = Pxy_z
//            //PressTensorGrad[n+12*Np] = Pyz_x
//            //PressTensorGrad[n+13*Np] = Pyz_y
//            //PressTensorGrad[n+14*Np] = Pyz_z
//            //PressTensorGrad[n+15*Np] = Pxz_x
//            //PressTensorGrad[n+16*Np] = Pxz_y
//            //PressTensorGrad[n+17*Np] = Pxz_z
//            Pxx_x = PressTensorGrad[0*Np+n];
//            Pyy_y = PressTensorGrad[4*Np+n];
//            Pzz_z = PressTensorGrad[8*Np+n];
//            Pxy_x = PressTensorGrad[9*Np+n];
//            Pxz_x = PressTensorGrad[15*Np+n];
//		    Pxy_y = PressTensorGrad[10*Np+n];
//		    Pyz_y = PressTensorGrad[13*Np+n];
//		    Pyz_z = PressTensorGrad[14*Np+n];
//		    Pxz_z = PressTensorGrad[17*Np+n];
//		    //............Compute the fluid-fluid force (gfx,gfy,gfz)...................................
//            //TODO double check if you need porosity as a fre-factor
//            Gff_x = porosity*px-(Pxx_x+Pxy_y+Pxz_z);
//            Gff_y = porosity*py-(Pxy_x+Pyy_y+Pyz_z);
//            Gff_z = porosity*pz-(Pxz_x+Pyz_y+Pzz_z);
//            // fluid-solid force
//            Gfs_x = (nA-nB)*SolidForce[n+0*Np];    
//            Gfs_y = (nA-nB)*SolidForce[n+1*Np];    
//            Gfs_z = (nA-nB)*SolidForce[n+2*Np];    
//
//			// local density
//			rho0=rhoA + 0.5*(1.0-phi)*(rhoB-rhoA);
//			// local relaxation time
//			tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);
//			rlx_setA = 1.f/tau;
//			rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
//			tau_eff=tauA_eff + 0.5*(1.0-phi)*(tauB_eff-tauA_eff);
//            mu_eff = (tau_eff-0.5)/3.f;//kinematic viscosity
//
//            //........................................................................
//            //					READ THE DISTRIBUTIONS
//            //		(read from opposite array due to previous swap operation)
//            //........................................................................
//            // q=0
//            fq = dist[n];
//            m1  = -30.0*fq;
//            m2  = 12.0*fq;
//
//            // q=1
//            nr1 = neighborList[n]; // neighbor 2 ( > 10Np => odd part of dist)
//            fq = dist[nr1]; // reading the f1 data into register fq
//            pressure = fq;
//            m1 -= 11.0*fq;
//            m2 -= 4.0*fq;
//            jx = fq;
//            m4 = -4.0*fq;
//            m9 = 2.0*fq;
//            m10 = -4.0*fq;
//
//            // q=2
//            nr2 = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
//            fq = dist[nr2];  // reading the f2 data into register fq
//            pressure += fq;
//            m1 -= 11.0*(fq);
//            m2 -= 4.0*(fq);
//            jx -= fq;
//            m4 += 4.0*(fq);
//            m9 += 2.0*(fq);
//            m10 -= 4.0*(fq);
//
//            // q=3
//            nr3 = neighborList[n+2*Np]; // neighbor 4
//            fq = dist[nr3];
//            pressure += fq;
//            m1 -= 11.0*fq;
//            m2 -= 4.0*fq;
//            jy = fq;
//            m6 = -4.0*fq;
//            m9 -= fq;
//            m10 += 2.0*fq;
//            m11 = fq;
//            m12 = -2.0*fq;
//
//            // q = 4
//            nr4 = neighborList[n+3*Np]; // neighbor 3
//            fq = dist[nr4];
//            pressure += fq;
//            m1 -= 11.0*fq;
//            m2 -= 4.0*fq;
//            jy -= fq;
//            m6 += 4.0*fq;
//            m9 -= fq;
//            m10 += 2.0*fq;
//            m11 += fq;
//            m12 -= 2.0*fq;
//
//            // q=5
//            nr5 = neighborList[n+4*Np];
//            fq = dist[nr5];
//            pressure += fq;
//            m1 -= 11.0*fq;
//            m2 -= 4.0*fq;
//            jz = fq;
//            m8 = -4.0*fq;
//            m9 -= fq;
//            m10 += 2.0*fq;
//            m11 -= fq;
//            m12 += 2.0*fq;
//
//            // q = 6
//            nr6 = neighborList[n+5*Np];
//            fq = dist[nr6];
//            pressure += fq;
//            m1 -= 11.0*fq;
//            m2 -= 4.0*fq;
//            jz -= fq;
//            m8 += 4.0*fq;
//            m9 -= fq;
//            m10 += 2.0*fq;
//            m11 -= fq;
//            m12 += 2.0*fq;
//
//            // q=7
//            nread = neighborList[n+6*Np];
//            fq = dist[nread];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jx += fq;
//            m4 += fq;
//            jy += fq;
//            m6 += fq;
//            m9  += fq;
//            m10 += fq;
//            m11 += fq;
//            m12 += fq;
//            m13 = fq;
//            m16 = fq;
//            m17 = -fq;
//
//            // q = 8
//            nread = neighborList[n+7*Np];
//            fq = dist[nread];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jx -= fq;
//            m4 -= fq;
//            jy -= fq;
//            m6 -= fq;
//            m9 += fq;
//            m10 += fq;
//            m11 += fq;
//            m12 += fq;
//            m13 += fq;
//            m16 -= fq;
//            m17 += fq;
//
//            // q=9
//            nread = neighborList[n+8*Np];
//            fq = dist[nread];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jx += fq;
//            m4 += fq;
//            jy -= fq;
//            m6 -= fq;
//            m9 += fq;
//            m10 += fq;
//            m11 += fq;
//            m12 += fq;
//            m13 -= fq;
//            m16 += fq;
//            m17 += fq;
//
//            // q = 10
//            nread = neighborList[n+9*Np];
//            fq = dist[nread];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jx -= fq;
//            m4 -= fq;
//            jy += fq;
//            m6 += fq;
//            m9 += fq;
//            m10 += fq;
//            m11 += fq;
//            m12 += fq;
//            m13 -= fq;
//            m16 -= fq;
//            m17 -= fq;
//
//            // q=11
//            nread = neighborList[n+10*Np];
//            fq = dist[nread];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jx += fq;
//            m4 += fq;
//            jz += fq;
//            m8 += fq;
//            m9 += fq;
//            m10 += fq;
//            m11 -= fq;
//            m12 -= fq;
//            m15 = fq;
//            m16 -= fq;
//            m18 = fq;
//
//            // q=12
//            nread = neighborList[n+11*Np];
//            fq = dist[nread];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jx -= fq;
//            m4 -= fq;
//            jz -= fq;
//            m8 -= fq;
//            m9 += fq;
//            m10 += fq;
//            m11 -= fq;
//            m12 -= fq;
//            m15 += fq;
//            m16 += fq;
//            m18 -= fq;
//
//            // q=13
//            nread = neighborList[n+12*Np];
//            fq = dist[nread];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jx += fq;
//            m4 += fq;
//            jz -= fq;
//            m8 -= fq;
//            m9 += fq;
//            m10 += fq;
//            m11 -= fq;
//            m12 -= fq;
//            m15 -= fq;
//            m16 -= fq;
//            m18 -= fq;
//
//            // q=14
//            nread = neighborList[n+13*Np];
//            fq = dist[nread];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jx -= fq;
//            m4 -= fq;
//            jz += fq;
//            m8 += fq;
//            m9 += fq;
//            m10 += fq;
//            m11 -= fq;
//            m12 -= fq;
//            m15 -= fq;
//            m16 += fq;
//            m18 += fq;
//
//            // q=15
//            nread = neighborList[n+14*Np];
//            fq = dist[nread];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jy += fq;
//            m6 += fq;
//            jz += fq;
//            m8 += fq;
//            m9 -= 2.0*fq;
//            m10 -= 2.0*fq;
//            m14 = fq;
//            m17 += fq;
//            m18 -= fq;
//
//            // q=16
//            nread = neighborList[n+15*Np];
//            fq = dist[nread];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jy -= fq;
//            m6 -= fq;
//            jz -= fq;
//            m8 -= fq;
//            m9 -= 2.0*fq;
//            m10 -= 2.0*fq;
//            m14 += fq;
//            m17 -= fq;
//            m18 += fq;
//
//            // q=17
//            nread = neighborList[n+16*Np];
//            fq = dist[nread];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jy += fq;
//            m6 += fq;
//            jz -= fq;
//            m8 -= fq;
//            m9 -= 2.0*fq;
//            m10 -= 2.0*fq;
//            m14 -= fq;
//            m17 += fq;
//            m18 += fq;
//
//            // q=18
//            nread = neighborList[n+17*Np];
//            fq = dist[nread];
//            pressure += fq;
//            m1 += 8.0*fq;
//            m2 += fq;
//            jy -= fq;
//            m6 -= fq;
//            jz += fq;
//            m8 += fq;
//            m9 -= 2.0*fq;
//            m10 -= 2.0*fq;
//            m14 -= fq;
//            m17 -= fq;
//            m18 -= fq;
//            //---------------------------------------------------------------------//
//
//            c0 = 0.5*(1.0+porosity*0.5*mu_eff/perm);
//            if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
//            //GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
//            c1 = porosity*0.5*GeoFun/sqrt(perm);
//            if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes
//
//            vx = jx/rho0+0.5*(porosity*Gx+Gff_x+Gfs_x);
//            vy = jy/rho0+0.5*(porosity*Gy+Gff_y+Gfs_y);
//            vz = jz/rho0+0.5*(porosity*Gz+Gff_z+Gfs_z);
//            v_mag=sqrt(vx*vx+vy*vy+vz*vz);
//            ux = vx/(c0+sqrt(c0*c0+c1*v_mag));
//            uy = vy/(c0+sqrt(c0*c0+c1*v_mag));
//            uz = vz/(c0+sqrt(c0*c0+c1*v_mag));
//            u_mag=sqrt(ux*ux+uy*uy+uz*uz);
//
//            //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
//            Fx = rho0*(-porosity*mu_eff/perm*ux - porosity*GeoFun/sqrt(perm)*u_mag*ux + porosity*Gx + Gff_x + Gfs_x);
//            Fy = rho0*(-porosity*mu_eff/perm*uy - porosity*GeoFun/sqrt(perm)*u_mag*uy + porosity*Gy + Gff_y + Gfs_y);
//            Fz = rho0*(-porosity*mu_eff/perm*uz - porosity*GeoFun/sqrt(perm)*u_mag*uz + porosity*Gz + Gff_z + Gfs_z);
//            if (porosity==1.0){
//                Fx=rho0*(Gx + Gff_x + Gfs_x);
//                Fy=rho0*(Gy + Gff_y + Gfs_y);
//                Fz=rho0*(Gz + Gff_z + Gfs_z);
//            }
//
//            //Calculate pressure for Incompressible-MRT model
//            pressure=0.5/porosity*(pressure-0.5*rho0*u_mag*u_mag/porosity);
//
////            //..............carry out relaxation process...............................................
////            m1 = m1 + rlx_setA*((-30*rho0+19*(ux*ux+uy*uy+uz*uz)/porosity + 57*pressure*porosity) - m1) 
////                    + (1-0.5*rlx_setA)*38*(Fx*ux+Fy*uy+Fz*uz)/porosity;
////            m2 = m2 + rlx_setA*((12*rho0 - 5.5*(ux*ux+uy*uy+uz*uz)/porosity-27*pressure*porosity) - m2)
////                    + (1-0.5*rlx_setA)*11*(-Fx*ux-Fy*uy-Fz*uz)/porosity;
////            jx = jx + Fx;
////            m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0) - m4)
////                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
////            jy = jy + Fy;
////            m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0) - m6)
////                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
////            jz = jz + Fz;
////            m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0) - m8)
////                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
////            m9 = m9 + rlx_setA*((rho0*(2*ux*ux-uy*uy-uz*uz)/porosity) - m9)
////                    + (1-0.5*rlx_setA)*(4*Fx*ux-2*Fy*uy-2*Fz*uz)/porosity;
////            m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10)
////                      + (1-0.5*rlx_setA)*(-2*Fx*ux+Fy*uy+Fz*uz)/porosity;
////            m11 = m11 + rlx_setA*((rho0*(uy*uy-uz*uz)/porosity) - m11)
////                      + (1-0.5*rlx_setA)*(2*Fy*uy-2*Fz*uz)/porosity;
////            m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12)
////                      + (1-0.5*rlx_setA)*(-Fy*uy+Fz*uz)/porosity;
////            m13 = m13 + rlx_setA*((rho0*ux*uy/porosity) - m13)
////                      + (1-0.5*rlx_setA)*(Fy*ux+Fx*uy)/porosity;
////            m14 = m14 + rlx_setA*((rho0*uy*uz/porosity) - m14)
////                      + (1-0.5*rlx_setA)*(Fz*uy+Fy*uz)/porosity;
////            m15 = m15 + rlx_setA*((rho0*ux*uz/porosity) - m15)
////                      + (1-0.5*rlx_setA)*(Fz*ux+Fx*uz)/porosity;
////            m16 = m16 + rlx_setB*( - m16);
////            m17 = m17 + rlx_setB*( - m17);
////            m18 = m18 + rlx_setB*( - m18);
////            //.......................................................................................................
//           
//            //-------------------- IMRT collison where body force has NO higher-order terms -------------//
//            //..............carry out relaxation process...............................................
//            m1 = m1 + rlx_setA*((-30*rho0+19*(ux*ux+uy*uy+uz*uz)/porosity + 57*pressure*porosity) - m1);
//            m2 = m2 + rlx_setA*((12*rho0 - 5.5*(ux*ux+uy*uy+uz*uz)/porosity-27*pressure*porosity) - m2);
//            jx = jx + Fx;
//            m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0) - m4)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
//            jy = jy + Fy;
//            m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0) - m6)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
//            jz = jz + Fz;
//            m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0) - m8)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
//            m9 = m9 + rlx_setA*((rho0*(2*ux*ux-uy*uy-uz*uz)/porosity) - m9);
//            m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10);
//            m11 = m11 + rlx_setA*((rho0*(uy*uy-uz*uz)/porosity) - m11);
//            m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12);
//            m13 = m13 + rlx_setA*((rho0*ux*uy/porosity) - m13);
//            m14 = m14 + rlx_setA*((rho0*uy*uz/porosity) - m14);
//            m15 = m15 + rlx_setA*((rho0*ux*uz/porosity) - m15);
//            m16 = m16 + rlx_setB*( - m16);
//            m17 = m17 + rlx_setB*( - m17);
//            m18 = m18 + rlx_setB*( - m18);
//            //.......................................................................................................
//
//
//            //.................inverse transformation......................................................
//            // q=0
//            fq = mrt_V1*rho0-mrt_V2*m1+mrt_V3*m2;
//            dist[n] = fq;
//
//            // q = 1
//            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10);
//            //nread = neighborList[n+Np];
//            dist[nr2] = fq;
//
//            // q=2
//            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10);
//            //nread = neighborList[n];
//            dist[nr1] = fq;
//
//            // q = 3
//            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
//            //nread = neighborList[n+3*Np];
//            dist[nr4] = fq;
//
//            // q = 4
//            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
//            //nread = neighborList[n+2*Np];
//            dist[nr3] = fq;
//
//            // q = 5
//            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
//            //nread = neighborList[n+5*Np];
//            dist[nr6] = fq;
//
//            // q = 6
//            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
//            //nread = neighborList[n+4*Np];
//            dist[nr5] = fq;
//
//            // q = 7
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m16-m17);
//            nread = neighborList[n+7*Np];
//            dist[nread] = fq;
//
//            // q = 8
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m17-m16);
//            nread = neighborList[n+6*Np];
//            dist[nread] = fq;
//
//            // q = 9
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13+0.125*(m16+m17);
//            nread = neighborList[n+9*Np];
//            dist[nread] = fq;
//
//            // q = 10
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13-0.125*(m16+m17);
//            nread = neighborList[n+8*Np];
//            dist[nread] = fq;
//
//            // q = 11
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m18-m16);
//            nread = neighborList[n+11*Np];
//            dist[nread] = fq;
//
//            // q = 12
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m16-m18);
//            nread = neighborList[n+10*Np];
//            dist[nread]= fq;
//
//            // q = 13
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12-0.25*m15-0.125*(m16+m18);
//            nread = neighborList[n+13*Np];
//            dist[nread] = fq;
//
//            // q= 14
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12-0.25*m15+0.125*(m16+m18);
//            nread = neighborList[n+12*Np];
//            dist[nread] = fq;
//
//            // q = 15
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18);
//            nread = neighborList[n+15*Np];
//            dist[nread] = fq;
//
//            // q = 16
//            fq =  mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17);
//            nread = neighborList[n+14*Np];
//            dist[nread] = fq;
//
//            // q = 17
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)-mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18);
//            nread = neighborList[n+17*Np];
//            dist[nread] = fq;
//
//            // q = 18
//            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)-mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18);
//            nread = neighborList[n+16*Np];
//            dist[nread] = fq;
//            //........................................................................
//
//            //Update velocity on device
//            Velocity[0*Np+n] = ux;
//            Velocity[1*Np+n] = uy;
//            Velocity[2*Np+n] = uz;
//            //Update pressure on device
//            Pressure[n] = pressure;
//
//            //-----------------------Mass transport------------------------//
//            // calcuale chemical potential
//            chem_a = lambdaA*(nA*nA*nA-1.5*nA*nA+0.5*nA)-0.25*kappaA*phi_lap;
//            chem_b = -lambdaB*(nB*nB*nB-1.5*nB*nB+0.5*nB)-0.25*kappaB*phi_lap;
//            rlx_massA = 3.f-sqrt(3.f);
//            rlx_massB = 3.f-sqrt(3.f);
//
//			//...............................................
//			// q = 0,2,4
//			// Cq = {1,0,0}, {0,1,0}, {0,0,1}
//			a1 = Aq[nr2];
//			b1 = Bq[nr2];
//			a2 = Aq[nr1];
//			b2 = Bq[nr1];
//			a1 = (1.0-rlx_massA)*a1+rlx_massA*(0.1111111111111111*4.5*(gamma*chem_a+nA*ux));
//			b1 = (1.0-rlx_massB)*b1+rlx_massB*(0.1111111111111111*4.5*(gamma*chem_b+nB*ux));
//			a2 = (1.0-rlx_massA)*a2+rlx_massA*(0.1111111111111111*4.5*(gamma*chem_a-nA*ux));
//			b2 = (1.0-rlx_massB)*b2+rlx_massB*(0.1111111111111111*4.5*(gamma*chem_b-nB*ux));
//
//			// q = 1
//			//nread = neighborList[n+Np];
//			Aq[nr2] = a1;
//			Bq[nr2] = b1;
//			// q=2
//			//nread = neighborList[n];
//			Aq[nr1] = a2;
//			Bq[nr1] = b2;
//
//			//...............................................
//			// Cq = {0,1,0}
//			a1 = Aq[nr4];
//			b1 = Bq[nr4];
//			a2 = Aq[nr3];
//			b2 = Bq[nr3];
//			a1 = (1.0-rlx_massA)*a1+rlx_massA*(0.1111111111111111*4.5*(gamma*chem_a+nA*uy));
//			b1 = (1.0-rlx_massB)*b1+rlx_massB*(0.1111111111111111*4.5*(gamma*chem_b+nB*uy));
//			a2 = (1.0-rlx_massA)*a2+rlx_massA*(0.1111111111111111*4.5*(gamma*chem_a-nA*uy));
//			b2 = (1.0-rlx_massB)*b2+rlx_massB*(0.1111111111111111*4.5*(gamma*chem_b-nB*uy));
//
//			// q = 3
//			//nread = neighborList[n+3*Np];
//			Aq[nr4] = a1;
//			Bq[nr4] = b1;
//			// q = 4
//			//nread = neighborList[n+2*Np];
//			Aq[nr3] = a2;
//			Bq[nr3] = b2;
//
//			//...............................................
//			// q = 4
//			// Cq = {0,0,1}
//			a1 = Aq[nr6];
//			b1 = Bq[nr6];
//			a2 = Aq[nr5];
//			b2 = Bq[nr5];
//			a1 = (1.0-rlx_massA)*a1+rlx_massA*(0.1111111111111111*4.5*(gamma*chem_a+nA*uz));
//			b1 = (1.0-rlx_massB)*b1+rlx_massB*(0.1111111111111111*4.5*(gamma*chem_b+nB*uz));
//			a2 = (1.0-rlx_massA)*a2+rlx_massA*(0.1111111111111111*4.5*(gamma*chem_a-nA*uz));
//			b2 = (1.0-rlx_massB)*b2+rlx_massB*(0.1111111111111111*4.5*(gamma*chem_b-nB*uz));
//
//			// q = 5
//			//nread = neighborList[n+5*Np];
//			Aq[nr6] = a1;
//			Bq[nr6] = b1;
//			// q = 6
//			//nread = neighborList[n+4*Np];
//			Aq[nr5] = a2;
//			Bq[nr5] = b2;
//			//...............................................
//
//			// Instantiate mass transport distributions
//			// Stationary value - distribution 0
//            a1=Aq[n];
//            b1=Bq[n];
//			Aq[n] = (1.0-rlx_massA)*a1+rlx_massA*(nA-3.0*gamma*chem_a);
//			Bq[n] = (1.0-rlx_massB)*b1+rlx_massB*(nB-3.0*gamma*chem_b);
//
//
//		}
//	}
//}

__global__ void dvc_ScaLBL_D3Q19_AAodd_GreyscaleColorChem(int *neighborList, double *dist, double *Cq, double *Phi, double *Den,double *SolidForce, int start, int finish, int Np,
                double tauA,double tauB,double tauA_eff,double tauB_eff,double rhoA,double rhoB,double gamma,double kappaA,double kappaB,double lambdaA,double lambdaB,
                double Gx, double Gy, double Gz,
                double *Poros,double *Perm, double *Velocity,double *Pressure,double *PressureGrad,double *PressTensorGrad,double *PhiLap){

	int n, nread, nr1,nr2,nr3,nr4,nr5,nr6;
	double vx,vy,vz,v_mag;
    double ux,uy,uz,u_mag;
    double pressure;//defined for this incompressible model
	// conserved momemnts
	double jx,jy,jz;
	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
    double fq;
    // currently disable 'GeoFun'
    double GeoFun=0.0;//geometric function from Guo's PRE 66, 036304 (2002)
    double porosity;
    double perm;//voxel permeability
    double c0, c1; //Guo's model parameters
    double Fx, Fy, Fz;//The total body force including Brinkman force and user-specified (Gx,Gy,Gz)
	double tau,tau_eff,rlx_setA,rlx_setB;
    double mu_eff;//effective kinematic viscosity for Darcy term
    double rho0;
    double phi;
    double phi_lap;//laplacian of phase field
    double nA,nB;
	//double a1,b1,a2,b2;
    double Gfs_x,Gfs_y,Gfs_z;
    double Gff_x,Gff_y,Gff_z;
    double chem;
    //double rlx_massA,rlx_massB;
    double rlx_phi;
    double a1,a2;//PDF of phase field
    // *---------------------------------Pressure Tensor Gradient------------------------------------*//
    double Pxx_x,Pyy_y,Pzz_z;
    double Pxy_x,Pxy_y;
    double Pyz_y,Pyz_z;
    double Pxz_x,Pxz_z;
    double px,py,pz; //pressure gradient

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

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
	    //........Get 1-D index for this thread....................
	    n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

		if ( n<finish ){		

			// read the component number densities
			nA = Den[n];
			nB = Den[Np + n];
            // read phase field
            phi = Phi[n];
            // load laplacian of phase field
            phi_lap = PhiLap[n];
            // Load voxel porosity and perm
            porosity = Poros[n];
            // use local saturation as an estimation of effective relperm values
            perm = Perm[n]*nA/(nA+nB)*int(phi>0.0)+Perm[n]*nB/(nA+nB)*int(phi<0.0);

            //Load pressure gradient
            px=PressureGrad[0*Np+n];
            py=PressureGrad[1*Np+n];
            pz=PressureGrad[2*Np+n];

            //Load pressure tensor gradient
            //For reference full list of PressTensorGrad
            //PressTensorGrad[n+0*Np]  = Pxx_x
            //PressTensorGrad[n+1*Np]  = Pxx_y
            //PressTensorGrad[n+2*Np]  = Pxx_z
            //PressTensorGrad[n+3*Np]  = Pyy_x
            //PressTensorGrad[n+4*Np]  = Pyy_y
            //PressTensorGrad[n+5*Np]  = Pyy_z
            //PressTensorGrad[n+6*Np]  = Pzz_x
            //PressTensorGrad[n+7*Np]  = Pzz_y
            //PressTensorGrad[n+8*Np]  = Pzz_z
            //PressTensorGrad[n+9*Np]  = Pxy_x
            //PressTensorGrad[n+10*Np] = Pxy_y
            //PressTensorGrad[n+11*Np] = Pxy_z
            //PressTensorGrad[n+12*Np] = Pyz_x
            //PressTensorGrad[n+13*Np] = Pyz_y
            //PressTensorGrad[n+14*Np] = Pyz_z
            //PressTensorGrad[n+15*Np] = Pxz_x
            //PressTensorGrad[n+16*Np] = Pxz_y
            //PressTensorGrad[n+17*Np] = Pxz_z
            Pxx_x = PressTensorGrad[0*Np+n];
            Pyy_y = PressTensorGrad[4*Np+n];
            Pzz_z = PressTensorGrad[8*Np+n];
            Pxy_x = PressTensorGrad[9*Np+n];
            Pxz_x = PressTensorGrad[15*Np+n];
		    Pxy_y = PressTensorGrad[10*Np+n];
		    Pyz_y = PressTensorGrad[13*Np+n];
		    Pyz_z = PressTensorGrad[14*Np+n];
		    Pxz_z = PressTensorGrad[17*Np+n];
		    //............Compute the fluid-fluid force (gfx,gfy,gfz)...................................
            //TODO double check if you need porosity as a fre-factor
            Gff_x = porosity*px-(Pxx_x+Pxy_y+Pxz_z);
            Gff_y = porosity*py-(Pxy_x+Pyy_y+Pyz_z);
            Gff_z = porosity*pz-(Pxz_x+Pyz_y+Pzz_z);
            // fluid-solid force
            Gfs_x = (nA-nB)*SolidForce[n+0*Np];    
            Gfs_y = (nA-nB)*SolidForce[n+1*Np];    
            Gfs_z = (nA-nB)*SolidForce[n+2*Np];    

			// local density
			rho0=rhoA + 0.5*(1.0-phi)*(rhoB-rhoA);
			// local relaxation time
			tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);
			rlx_setA = 1.f/tau;
			rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
			tau_eff=tauA_eff + 0.5*(1.0-phi)*(tauB_eff-tauA_eff);
            mu_eff = (tau_eff-0.5)/3.f;//kinematic viscosity

            //........................................................................
            //					READ THE DISTRIBUTIONS
            //		(read from opposite array due to previous swap operation)
            //........................................................................
            // q=0
            fq = dist[n];
            m1  = -30.0*fq;
            m2  = 12.0*fq;

            // q=1
            nr1 = neighborList[n]; // neighbor 2 ( > 10Np => odd part of dist)
            fq = dist[nr1]; // reading the f1 data into register fq
            pressure = fq;
            m1 -= 11.0*fq;
            m2 -= 4.0*fq;
            jx = fq;
            m4 = -4.0*fq;
            m9 = 2.0*fq;
            m10 = -4.0*fq;

            // q=2
            nr2 = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
            fq = dist[nr2];  // reading the f2 data into register fq
            pressure += fq;
            m1 -= 11.0*(fq);
            m2 -= 4.0*(fq);
            jx -= fq;
            m4 += 4.0*(fq);
            m9 += 2.0*(fq);
            m10 -= 4.0*(fq);

            // q=3
            nr3 = neighborList[n+2*Np]; // neighbor 4
            fq = dist[nr3];
            pressure += fq;
            m1 -= 11.0*fq;
            m2 -= 4.0*fq;
            jy = fq;
            m6 = -4.0*fq;
            m9 -= fq;
            m10 += 2.0*fq;
            m11 = fq;
            m12 = -2.0*fq;

            // q = 4
            nr4 = neighborList[n+3*Np]; // neighbor 3
            fq = dist[nr4];
            pressure += fq;
            m1 -= 11.0*fq;
            m2 -= 4.0*fq;
            jy -= fq;
            m6 += 4.0*fq;
            m9 -= fq;
            m10 += 2.0*fq;
            m11 += fq;
            m12 -= 2.0*fq;

            // q=5
            nr5 = neighborList[n+4*Np];
            fq = dist[nr5];
            pressure += fq;
            m1 -= 11.0*fq;
            m2 -= 4.0*fq;
            jz = fq;
            m8 = -4.0*fq;
            m9 -= fq;
            m10 += 2.0*fq;
            m11 -= fq;
            m12 += 2.0*fq;

            // q = 6
            nr6 = neighborList[n+5*Np];
            fq = dist[nr6];
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            nread = neighborList[n+16*Np];
            fq = dist[nread];
            pressure += fq;
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
            pressure += fq;
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
            //---------------------------------------------------------------------//

            c0 = 0.5*(1.0+porosity*0.5*mu_eff/perm);
            if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
            //GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
            c1 = porosity*0.5*GeoFun/sqrt(perm);
            if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes

            vx = jx/rho0+0.5*(porosity*Gx+Gff_x+Gfs_x);
            vy = jy/rho0+0.5*(porosity*Gy+Gff_y+Gfs_y);
            vz = jz/rho0+0.5*(porosity*Gz+Gff_z+Gfs_z);
            v_mag=sqrt(vx*vx+vy*vy+vz*vz);
            ux = vx/(c0+sqrt(c0*c0+c1*v_mag));
            uy = vy/(c0+sqrt(c0*c0+c1*v_mag));
            uz = vz/(c0+sqrt(c0*c0+c1*v_mag));
            u_mag=sqrt(ux*ux+uy*uy+uz*uz);

            //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
            Fx = rho0*(-porosity*mu_eff/perm*ux - porosity*GeoFun/sqrt(perm)*u_mag*ux + porosity*Gx + Gff_x + Gfs_x);
            Fy = rho0*(-porosity*mu_eff/perm*uy - porosity*GeoFun/sqrt(perm)*u_mag*uy + porosity*Gy + Gff_y + Gfs_y);
            Fz = rho0*(-porosity*mu_eff/perm*uz - porosity*GeoFun/sqrt(perm)*u_mag*uz + porosity*Gz + Gff_z + Gfs_z);
            if (porosity==1.0){
                Fx=rho0*(Gx + Gff_x + Gfs_x);
                Fy=rho0*(Gy + Gff_y + Gfs_y);
                Fz=rho0*(Gz + Gff_z + Gfs_z);
            }

            //Calculate pressure for Incompressible-MRT model
            pressure=0.5/porosity*(pressure-0.5*rho0*u_mag*u_mag/porosity);

//            //..............carry out relaxation process...............................................
//            m1 = m1 + rlx_setA*((-30*rho0+19*(ux*ux+uy*uy+uz*uz)/porosity + 57*pressure*porosity) - m1) 
//                    + (1-0.5*rlx_setA)*38*(Fx*ux+Fy*uy+Fz*uz)/porosity;
//            m2 = m2 + rlx_setA*((12*rho0 - 5.5*(ux*ux+uy*uy+uz*uz)/porosity-27*pressure*porosity) - m2)
//                    + (1-0.5*rlx_setA)*11*(-Fx*ux-Fy*uy-Fz*uz)/porosity;
//            jx = jx + Fx;
//            m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0) - m4)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
//            jy = jy + Fy;
//            m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0) - m6)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
//            jz = jz + Fz;
//            m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0) - m8)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
//            m9 = m9 + rlx_setA*((rho0*(2*ux*ux-uy*uy-uz*uz)/porosity) - m9)
//                    + (1-0.5*rlx_setA)*(4*Fx*ux-2*Fy*uy-2*Fz*uz)/porosity;
//            m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10)
//                      + (1-0.5*rlx_setA)*(-2*Fx*ux+Fy*uy+Fz*uz)/porosity;
//            m11 = m11 + rlx_setA*((rho0*(uy*uy-uz*uz)/porosity) - m11)
//                      + (1-0.5*rlx_setA)*(2*Fy*uy-2*Fz*uz)/porosity;
//            m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12)
//                      + (1-0.5*rlx_setA)*(-Fy*uy+Fz*uz)/porosity;
//            m13 = m13 + rlx_setA*((rho0*ux*uy/porosity) - m13)
//                      + (1-0.5*rlx_setA)*(Fy*ux+Fx*uy)/porosity;
//            m14 = m14 + rlx_setA*((rho0*uy*uz/porosity) - m14)
//                      + (1-0.5*rlx_setA)*(Fz*uy+Fy*uz)/porosity;
//            m15 = m15 + rlx_setA*((rho0*ux*uz/porosity) - m15)
//                      + (1-0.5*rlx_setA)*(Fz*ux+Fx*uz)/porosity;
//            m16 = m16 + rlx_setB*( - m16);
//            m17 = m17 + rlx_setB*( - m17);
//            m18 = m18 + rlx_setB*( - m18);
//            //.......................................................................................................
           
            //-------------------- IMRT collison where body force has NO higher-order terms -------------//
            //..............carry out relaxation process...............................................
            m1 = m1 + rlx_setA*((-30*rho0+19*(ux*ux+uy*uy+uz*uz)/porosity + 57*pressure*porosity) - m1);
            m2 = m2 + rlx_setA*((12*rho0 - 5.5*(ux*ux+uy*uy+uz*uz)/porosity-27*pressure*porosity) - m2);
            jx = jx + Fx;
            m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0) - m4)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
            jy = jy + Fy;
            m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0) - m6)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
            jz = jz + Fz;
            m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0) - m8)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
            m9 = m9 + rlx_setA*((rho0*(2*ux*ux-uy*uy-uz*uz)/porosity) - m9);
            m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10);
            m11 = m11 + rlx_setA*((rho0*(uy*uy-uz*uz)/porosity) - m11);
            m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12);
            m13 = m13 + rlx_setA*((rho0*ux*uy/porosity) - m13);
            m14 = m14 + rlx_setA*((rho0*uy*uz/porosity) - m14);
            m15 = m15 + rlx_setA*((rho0*ux*uz/porosity) - m15);
            m16 = m16 + rlx_setB*( - m16);
            m17 = m17 + rlx_setB*( - m17);
            m18 = m18 + rlx_setB*( - m18);
            //.......................................................................................................


            //.................inverse transformation......................................................
            // q=0
            fq = mrt_V1*rho0-mrt_V2*m1+mrt_V3*m2;
            dist[n] = fq;

            // q = 1
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10);
            //nread = neighborList[n+Np];
            dist[nr2] = fq;

            // q=2
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10);
            //nread = neighborList[n];
            dist[nr1] = fq;

            // q = 3
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
            //nread = neighborList[n+3*Np];
            dist[nr4] = fq;

            // q = 4
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
            //nread = neighborList[n+2*Np];
            dist[nr3] = fq;

            // q = 5
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
            //nread = neighborList[n+5*Np];
            dist[nr6] = fq;

            // q = 6
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
            //nread = neighborList[n+4*Np];
            dist[nr5] = fq;

            // q = 7
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m16-m17);
            nread = neighborList[n+7*Np];
            dist[nread] = fq;

            // q = 8
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m17-m16);
            nread = neighborList[n+6*Np];
            dist[nread] = fq;

            // q = 9
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13+0.125*(m16+m17);
            nread = neighborList[n+9*Np];
            dist[nread] = fq;

            // q = 10
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13-0.125*(m16+m17);
            nread = neighborList[n+8*Np];
            dist[nread] = fq;

            // q = 11
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m18-m16);
            nread = neighborList[n+11*Np];
            dist[nread] = fq;

            // q = 12
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m16-m18);
            nread = neighborList[n+10*Np];
            dist[nread]= fq;

            // q = 13
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12-0.25*m15-0.125*(m16+m18);
            nread = neighborList[n+13*Np];
            dist[nread] = fq;

            // q= 14
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12-0.25*m15+0.125*(m16+m18);
            nread = neighborList[n+12*Np];
            dist[nread] = fq;

            // q = 15
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18);
            nread = neighborList[n+15*Np];
            dist[nread] = fq;

            // q = 16
            fq =  mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17);
            nread = neighborList[n+14*Np];
            dist[nread] = fq;

            // q = 17
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)-mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18);
            nread = neighborList[n+17*Np];
            dist[nread] = fq;

            // q = 18
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)-mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18);
            nread = neighborList[n+16*Np];
            dist[nread] = fq;
            //........................................................................

            //Update velocity on device
            Velocity[0*Np+n] = ux;
            Velocity[1*Np+n] = uy;
            Velocity[2*Np+n] = uz;
            //Update pressure on device
            Pressure[n] = pressure;

            //-----------------------Mass transport------------------------//
            // calcuale chemical potential
            chem = lambdaA*(nA*nA*nA-1.5*nA*nA+0.5*nA)-lambdaB*(nB*nB*nB-1.5*nB*nB+0.5*nB)-0.25*(kappaA+kappaB)*phi_lap;
            //rlx_phi = 3.f-sqrt(3.f);
            rlx_phi = 1.0;

			//...............................................
			// q = 0,2,4
			// Cq = {1,0,0}, {0,1,0}, {0,0,1}
			a1 = Cq[nr2];
			a2 = Cq[nr1];
			a1 = (1.0-rlx_phi)*a1+rlx_phi*(0.1111111111111111*4.5*(gamma*chem+phi*ux));
			a2 = (1.0-rlx_phi)*a2+rlx_phi*(0.1111111111111111*4.5*(gamma*chem-phi*ux));

			// q = 1
			//nread = neighborList[n+Np];
			Cq[nr2] = a1;
			// q=2
			//nread = neighborList[n];
			Cq[nr1] = a2;

			//...............................................
			// Cq = {0,1,0}
			a1 = Cq[nr4];
			a2 = Cq[nr3];
			a1 = (1.0-rlx_phi)*a1+rlx_phi*(0.1111111111111111*4.5*(gamma*chem+phi*uy));
			a2 = (1.0-rlx_phi)*a2+rlx_phi*(0.1111111111111111*4.5*(gamma*chem-phi*uy));

			// q = 3
			//nread = neighborList[n+3*Np];
			Cq[nr4] = a1;
			// q = 4
			//nread = neighborList[n+2*Np];
			Cq[nr3] = a2;

			//...............................................
			// q = 4
			// Cq = {0,0,1}
			a1 = Cq[nr6];
			a2 = Cq[nr5];
			a1 = (1.0-rlx_phi)*a1+rlx_phi*(0.1111111111111111*4.5*(gamma*chem+phi*uz));
			a2 = (1.0-rlx_phi)*a2+rlx_phi*(0.1111111111111111*4.5*(gamma*chem-phi*uz));

			// q = 5
			//nread = neighborList[n+5*Np];
			Cq[nr6] = a1;
			// q = 6
			//nread = neighborList[n+4*Np];
			Cq[nr5] = a2;
			//...............................................

			// Instantiate mass transport distributions
			// Stationary value - distribution 0
            a1=Cq[n];
			Cq[n] = (1.0-rlx_phi)*a1+rlx_phi*(a1-3.0*gamma*chem);

		}
	}
}

__global__ void dvc_ScaLBL_D3Q19_AAeven_GreyscaleColorChem(double *dist, double *Cq, double *Phi, double *Den,double *SolidForce, int start, int finish, int Np,
                double tauA,double tauB,double tauA_eff,double tauB_eff,double rhoA,double rhoB,double gamma,double kappaA,double kappaB,double lambdaA,double lambdaB,
                double Gx, double Gy, double Gz,
                double *Poros,double *Perm, double *Velocity,double *Pressure,double *PressureGrad,double *PressTensorGrad,double *PhiLap){
	int n;
	double vx,vy,vz,v_mag;
    double ux,uy,uz,u_mag;
    double pressure;//defined for this incompressible model
	// conserved momemnts
	double jx,jy,jz;
	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
    double fq;
    // currently disable 'GeoFun'
    double GeoFun=0.0;//geometric function from Guo's PRE 66, 036304 (2002)
    double porosity;
    double perm;//voxel permeability
    double c0, c1; //Guo's model parameters
    double Fx, Fy, Fz;//The total body force including Brinkman force and user-specified (Gx,Gy,Gz)
	double tau,tau_eff,rlx_setA,rlx_setB;
    double mu_eff;//effective kinematic viscosity for Darcy term
    double rho0;
    double phi;
    double phi_lap;//laplacian of phase field
    double nA,nB;
	//double a1,b1,a2,b2;
    double Gfs_x,Gfs_y,Gfs_z;
    double Gff_x,Gff_y,Gff_z;
    double chem;
    //double rlx_massA,rlx_massB;
    double rlx_phi;
    double a1,a2;//PDF of phase field
    // *---------------------------------Pressure Tensor Gradient------------------------------------*//
    double Pxx_x,Pyy_y,Pzz_z;
    double Pxy_x,Pxy_y;
    double Pyz_y,Pyz_z;
    double Pxz_x,Pxz_z;
    double px,py,pz; //pressure gradient


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


	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
	    //........Get 1-D index for this thread....................
	    n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

		if ( n<finish ){

			// read the component number densities
            // TODO you can eliminate this, get nA and nB from phi
			nA = Den[n];
			nB = Den[Np + n];
            // read phase field
            phi = Phi[n];
            // load laplacian of phase field
            phi_lap = PhiLap[n];
            // Load voxel porosity and perm
            porosity = Poros[n];
            // use local saturation as an estimation of effective relperm values
            perm = Perm[n]*nA/(nA+nB)*int(phi>0.0)+Perm[n]*nB/(nA+nB)*int(phi<0.0);

            //Load pressure gradient
            px=PressureGrad[0*Np+n];
            py=PressureGrad[1*Np+n];
            pz=PressureGrad[2*Np+n];

            //Load pressure tensor gradient
            //For reference full list of PressTensorGrad
            //PressTensorGrad[n+0*Np]  = Pxx_x
            //PressTensorGrad[n+1*Np]  = Pxx_y
            //PressTensorGrad[n+2*Np]  = Pxx_z
            //PressTensorGrad[n+3*Np]  = Pyy_x
            //PressTensorGrad[n+4*Np]  = Pyy_y
            //PressTensorGrad[n+5*Np]  = Pyy_z
            //PressTensorGrad[n+6*Np]  = Pzz_x
            //PressTensorGrad[n+7*Np]  = Pzz_y
            //PressTensorGrad[n+8*Np]  = Pzz_z
            //PressTensorGrad[n+9*Np]  = Pxy_x
            //PressTensorGrad[n+10*Np] = Pxy_y
            //PressTensorGrad[n+11*Np] = Pxy_z
            //PressTensorGrad[n+12*Np] = Pyz_x
            //PressTensorGrad[n+13*Np] = Pyz_y
            //PressTensorGrad[n+14*Np] = Pyz_z
            //PressTensorGrad[n+15*Np] = Pxz_x
            //PressTensorGrad[n+16*Np] = Pxz_y
            //PressTensorGrad[n+17*Np] = Pxz_z
            Pxx_x = PressTensorGrad[0*Np+n];
            Pyy_y = PressTensorGrad[4*Np+n];
            Pzz_z = PressTensorGrad[8*Np+n];
            Pxy_x = PressTensorGrad[9*Np+n];
            Pxz_x = PressTensorGrad[15*Np+n];
		    Pxy_y = PressTensorGrad[10*Np+n];
		    Pyz_y = PressTensorGrad[13*Np+n];
		    Pyz_z = PressTensorGrad[14*Np+n];
		    Pxz_z = PressTensorGrad[17*Np+n];
		    //............Compute the fluid-fluid force (gfx,gfy,gfz)...................................
            //TODO double check if you need porosity as a fre-factor
            Gff_x = porosity*px-(Pxx_x+Pxy_y+Pxz_z);
            Gff_y = porosity*py-(Pxy_x+Pyy_y+Pyz_z);
            Gff_z = porosity*pz-(Pxz_x+Pyz_y+Pzz_z);
            // fluid-solid force
            Gfs_x = (nA-nB)*SolidForce[n+0*Np];    
            Gfs_y = (nA-nB)*SolidForce[n+1*Np];    
            Gfs_z = (nA-nB)*SolidForce[n+2*Np];    

			// local density
			rho0=rhoA + 0.5*(1.0-phi)*(rhoB-rhoA);
			// local relaxation time
			tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);
			rlx_setA = 1.f/tau;
			rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
			tau_eff=tauA_eff + 0.5*(1.0-phi)*(tauB_eff-tauA_eff);
            mu_eff = (tau_eff-0.5)/3.f;//kinematic viscosity


            //........................................................................
            //					READ THE DISTRIBUTIONS
            //		(read from opposite array due to previous swap operation)
            //........................................................................
            // q=0
            fq = dist[n];
            m1  = -30.0*fq;
            m2  = 12.0*fq;

            // q=1
            fq = dist[2*Np+n];
            pressure = fq;
            m1 -= 11.0*fq;
            m2 -= 4.0*fq;
            jx = fq;
            m4 = -4.0*fq;
            m9 = 2.0*fq;
            m10 = -4.0*fq;

            // f2 = dist[10*Np+n];
            fq = dist[1*Np+n];
            pressure += fq;
            m1 -= 11.0*(fq);
            m2 -= 4.0*(fq);
            jx -= fq;
            m4 += 4.0*(fq);
            m9 += 2.0*(fq);
            m10 -= 4.0*(fq);

            // q=3
            fq = dist[4*Np+n];
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            pressure += fq;
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
            //---------------------------------------------------------------------//

            c0 = 0.5*(1.0+porosity*0.5*mu_eff/perm);
            if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
            //GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
            c1 = porosity*0.5*GeoFun/sqrt(perm);
            if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes

            vx = jx/rho0+0.5*(porosity*Gx+Gff_x+Gfs_x);
            vy = jy/rho0+0.5*(porosity*Gy+Gff_y+Gfs_y);
            vz = jz/rho0+0.5*(porosity*Gz+Gff_z+Gfs_z);
            v_mag=sqrt(vx*vx+vy*vy+vz*vz);
            ux = vx/(c0+sqrt(c0*c0+c1*v_mag));
            uy = vy/(c0+sqrt(c0*c0+c1*v_mag));
            uz = vz/(c0+sqrt(c0*c0+c1*v_mag));
            u_mag=sqrt(ux*ux+uy*uy+uz*uz);

            //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
            Fx = rho0*(-porosity*mu_eff/perm*ux - porosity*GeoFun/sqrt(perm)*u_mag*ux + porosity*Gx + Gff_x + Gfs_x);
            Fy = rho0*(-porosity*mu_eff/perm*uy - porosity*GeoFun/sqrt(perm)*u_mag*uy + porosity*Gy + Gff_y + Gfs_y);
            Fz = rho0*(-porosity*mu_eff/perm*uz - porosity*GeoFun/sqrt(perm)*u_mag*uz + porosity*Gz + Gff_z + Gfs_z);
            if (porosity==1.0){
                Fx=rho0*(Gx + Gff_x + Gfs_x);
                Fy=rho0*(Gy + Gff_y + Gfs_y);
                Fz=rho0*(Gz + Gff_z + Gfs_z);
            }

            //Calculate pressure for Incompressible-MRT model
            pressure=0.5/porosity*(pressure-0.5*rho0*u_mag*u_mag/porosity);

//            //..............carry out relaxation process...............................................
//            m1 = m1 + rlx_setA*((-30*rho0+19*(ux*ux+uy*uy+uz*uz)/porosity + 57*pressure*porosity) - m1) 
//                    + (1-0.5*rlx_setA)*38*(Fx*ux+Fy*uy+Fz*uz)/porosity;
//            m2 = m2 + rlx_setA*((12*rho0 - 5.5*(ux*ux+uy*uy+uz*uz)/porosity-27*pressure*porosity) - m2)
//                    + (1-0.5*rlx_setA)*11*(-Fx*ux-Fy*uy-Fz*uz)/porosity;
//            jx = jx + Fx;
//            m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0) - m4)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
//            jy = jy + Fy;
//            m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0) - m6)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
//            jz = jz + Fz;
//            m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0) - m8)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
//            m9 = m9 + rlx_setA*((rho0*(2*ux*ux-uy*uy-uz*uz)/porosity) - m9)
//                    + (1-0.5*rlx_setA)*(4*Fx*ux-2*Fy*uy-2*Fz*uz)/porosity;
//            m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10)
//                      + (1-0.5*rlx_setA)*(-2*Fx*ux+Fy*uy+Fz*uz)/porosity;
//            m11 = m11 + rlx_setA*((rho0*(uy*uy-uz*uz)/porosity) - m11)
//                      + (1-0.5*rlx_setA)*(2*Fy*uy-2*Fz*uz)/porosity;
//            m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12)
//                      + (1-0.5*rlx_setA)*(-Fy*uy+Fz*uz)/porosity;
//            m13 = m13 + rlx_setA*((rho0*ux*uy/porosity) - m13)
//                      + (1-0.5*rlx_setA)*(Fy*ux+Fx*uy)/porosity;
//            m14 = m14 + rlx_setA*((rho0*uy*uz/porosity) - m14)
//                      + (1-0.5*rlx_setA)*(Fz*uy+Fy*uz)/porosity;
//            m15 = m15 + rlx_setA*((rho0*ux*uz/porosity) - m15)
//                      + (1-0.5*rlx_setA)*(Fz*ux+Fx*uz)/porosity;
//            m16 = m16 + rlx_setB*( - m16);
//            m17 = m17 + rlx_setB*( - m17);
//            m18 = m18 + rlx_setB*( - m18);
//            //.......................................................................................................

            //-------------------- IMRT collison where body force has NO higher-order terms -------------//
            //..............carry out relaxation process...............................................
            m1 = m1 + rlx_setA*((-30*rho0+19*(ux*ux+uy*uy+uz*uz)/porosity + 57*pressure*porosity) - m1);
            m2 = m2 + rlx_setA*((12*rho0 - 5.5*(ux*ux+uy*uy+uz*uz)/porosity-27*pressure*porosity) - m2);
            jx = jx + Fx;
            m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0) - m4)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
            jy = jy + Fy;
            m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0) - m6)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
            jz = jz + Fz;
            m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0) - m8)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
            m9 = m9 + rlx_setA*((rho0*(2*ux*ux-uy*uy-uz*uz)/porosity) - m9);
            m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10);
            m11 = m11 + rlx_setA*((rho0*(uy*uy-uz*uz)/porosity) - m11);
            m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12);
            m13 = m13 + rlx_setA*((rho0*ux*uy/porosity) - m13);
            m14 = m14 + rlx_setA*((rho0*uy*uz/porosity) - m14);
            m15 = m15 + rlx_setA*((rho0*ux*uz/porosity) - m15);
            m16 = m16 + rlx_setB*( - m16);
            m17 = m17 + rlx_setB*( - m17);
            m18 = m18 + rlx_setB*( - m18);
            //.......................................................................................................

            //.................inverse transformation......................................................
            // q=0
            fq = mrt_V1*rho0-mrt_V2*m1+mrt_V3*m2;
            dist[n] = fq;

            // q = 1
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10);
            dist[1*Np+n] = fq;

            // q=2
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10);
            dist[2*Np+n] = fq;

            // q = 3
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
            dist[3*Np+n] = fq;

            // q = 4
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
            dist[4*Np+n] = fq;

            // q = 5
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
            dist[5*Np+n] = fq;

            // q = 6
            fq = mrt_V1*rho0-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
            dist[6*Np+n] = fq;

            // q = 7
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m16-m17);
            dist[7*Np+n] = fq;

            // q = 8
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m17-m16);
            dist[8*Np+n] = fq;

            // q = 9
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13+0.125*(m16+m17);
            dist[9*Np+n] = fq;

            // q = 10
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13-0.125*(m16+m17);
            dist[10*Np+n] = fq;

            // q = 11
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m18-m16);
            dist[11*Np+n] = fq;

            // q = 12
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m16-m18);
            dist[12*Np+n] = fq;

            // q = 13
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12-0.25*m15-0.125*(m16+m18);
            dist[13*Np+n] = fq;

            // q= 14
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12-0.25*m15+0.125*(m16+m18);
            dist[14*Np+n] = fq;

            // q = 15
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18);
            dist[15*Np+n] = fq;

            // q = 16
            fq =  mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17);
            dist[16*Np+n] = fq;

            // q = 17
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)-mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18);
            dist[17*Np+n] = fq;

            // q = 18
            fq = mrt_V1*rho0+mrt_V9*m1+mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)-mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18);
            dist[18*Np+n] = fq;
            //........................................................................

            //Update velocity on device
            Velocity[0*Np+n] = ux;
            Velocity[1*Np+n] = uy;
            Velocity[2*Np+n] = uz;
            //Update pressure on device
            Pressure[n] = pressure;

            //-----------------------Mass transport------------------------//
            // calcuale chemical potential
            chem = lambdaA*(nA*nA*nA-1.5*nA*nA+0.5*nA)-lambdaB*(nB*nB*nB-1.5*nB*nB+0.5*nB)-0.25*(kappaA+kappaB)*phi_lap;
            //rlx_phi = 3.f-sqrt(3.f);
            rlx_phi = 1.0;

			//...............................................
			// q = 0,2,4
			// Cq = {1,0,0}, {0,1,0}, {0,0,1}
			a1 = Cq[1*Np+n];
			a2 = Cq[2*Np+n];
			a1 = (1.0-rlx_phi)*a1+rlx_phi*(0.1111111111111111*4.5*(gamma*chem+phi*ux));
			a2 = (1.0-rlx_phi)*a2+rlx_phi*(0.1111111111111111*4.5*(gamma*chem-phi*ux));

			Cq[1*Np+n] = a1;
			Cq[2*Np+n] = a2;

			//...............................................
			// q = 2
			// Cq = {0,1,0}
			a1 = Cq[3*Np+n];
			a2 = Cq[4*Np+n];
			a1 = (1.0-rlx_phi)*a1+rlx_phi*(0.1111111111111111*4.5*(gamma*chem+phi*uy));
			a2 = (1.0-rlx_phi)*a2+rlx_phi*(0.1111111111111111*4.5*(gamma*chem-phi*uy));

			Cq[3*Np+n] = a1;
			Cq[4*Np+n] = a2;
			//...............................................
			// q = 4
			// Cq = {0,0,1}
			a1 = Cq[5*Np+n];
			a2 = Cq[6*Np+n];
			a1 = (1.0-rlx_phi)*a1+rlx_phi*(0.1111111111111111*4.5*(gamma*chem+phi*uz));
			a2 = (1.0-rlx_phi)*a2+rlx_phi*(0.1111111111111111*4.5*(gamma*chem-phi*uz));

			Cq[5*Np+n] = a1;
			Cq[6*Np+n] = a2;
			//...............................................

			// Instantiate mass transport distributions
			// Stationary value - distribution 0
            a1=Cq[n];
			Cq[n] = (1.0-rlx_phi)*a1+rlx_phi*(a1-3.0*gamma*chem);
		}
	}
}

__global__ void dvc_ScaLBL_D3Q19_GreyColorIMRT_Init(double *dist, double *Den, double rhoA, double rhoB, int Np){
	int n;
	int S = Np/NBLOCKS/NTHREADS + 1;
    double phi;
    double nA,nB;
    double Den0;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		if (n<Np ){
            nA = Den[n];
            nB = Den[n+Np];
            phi = (nA-nB)/(nA+nB);
            Den0 = 0.5*(1.f+phi)*rhoA + 0.5*(1.f-phi)*rhoB;

			dist[n] = Den0 - 0.6666666666666667;
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
}

//__global__ void dvc_ScaLBL_D3Q7_GreyColorIMRT_Init(double *Den, double *Aq, double *Bq, double *Phi, int start, int finish, int Np){
//	int idx;
//    double nA,nB;
//
//	int S = Np/NBLOCKS/NTHREADS + 1;
//	for (int s=0; s<S; s++){
//		//........Get 1-D index for this thread....................
//		idx =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
//		if (idx<finish) {
//            nA = Den[idx];
//            nB = Den[idx+Np];
//
//			Aq[idx]=0.3333333333333333*nA;
//			Aq[Np+idx]=0.1111111111111111*nA;
//			Aq[2*Np+idx]=0.1111111111111111*nA;
//			Aq[3*Np+idx]=0.1111111111111111*nA;
//			Aq[4*Np+idx]=0.1111111111111111*nA;
//			Aq[5*Np+idx]=0.1111111111111111*nA;
//			Aq[6*Np+idx]=0.1111111111111111*nA;
//
//			Bq[idx]=0.3333333333333333*nB;
//			Bq[Np+idx]=0.1111111111111111*nB;
//			Bq[2*Np+idx]=0.1111111111111111*nB;
//			Bq[3*Np+idx]=0.1111111111111111*nB;
//			Bq[4*Np+idx]=0.1111111111111111*nB;
//			Bq[5*Np+idx]=0.1111111111111111*nB;
//			Bq[6*Np+idx]=0.1111111111111111*nB;
//
//            Phi[idx] = nA-nB;
//		}
//	}
//}

__global__ void dvc_ScaLBL_D3Q7_GreyColorIMRT_Init(double *Den, double *Cq, double *PhiLap, double gamma, double kappaA, double kappaB, double lambdaA, double lambdaB,
                int start, int finish, int Np){
	int idx;
    double nA,nB;
    double phi;
    double phi_lap;//laplacian of the phase field
    double chem;//chemical potential
	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		idx =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (idx<finish) {
            nA = Den[idx];
            nB = Den[idx+Np];
            phi = nA-nB;
            phi_lap = PhiLap[idx];
            chem = lambdaA*(nA*nA*nA-1.5*nA*nA+0.5*nA)-lambdaB*(nB*nB*nB-1.5*nB*nB+0.5*nB)-0.25*(kappaA+kappaB)*phi_lap;

			Cq[1*Np+idx]=0.5*gamma*chem;
			Cq[2*Np+idx]=0.5*gamma*chem;
			Cq[3*Np+idx]=0.5*gamma*chem;
			Cq[4*Np+idx]=0.5*gamma*chem;
			Cq[5*Np+idx]=0.5*gamma*chem;
			Cq[6*Np+idx]=0.5*gamma*chem;

			Cq[0*Np+idx]= phi - 3.0*gamma*chem;
		}
	}
}

__global__  void dvc_ScaLBL_D3Q7_AAodd_GreyscaleColorDensity(int *neighborList, double *Aq, double *Bq, double *Den, double *Phi, int start, int finish, int Np){
	int n,nread;
	double fq,nA,nB;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {
			//..........Compute the number density for each component ............
			// q=0
			fq = Aq[n];
			nA = fq;
			fq = Bq[n];
			nB = fq;
			
			// q=1
			nread = neighborList[n]; 
			fq = Aq[nread];
			nA += fq;
			fq = Bq[nread]; 
			nB += fq;
			
			// q=2
			nread = neighborList[n+Np]; 
			fq = Aq[nread];  
			nA += fq;
			fq = Bq[nread]; 
			nB += fq;
			
			// q=3
			nread = neighborList[n+2*Np]; 
			fq = Aq[nread];
			nA += fq;
			fq = Bq[nread];
			nB += fq;
			
			// q = 4
			nread = neighborList[n+3*Np]; 
			fq = Aq[nread];
			nA += fq;
			fq = Bq[nread];
			nB += fq;

			// q=5
			nread = neighborList[n+4*Np];
			fq = Aq[nread];
			nA += fq;
			fq = Bq[nread];
			nB += fq;
			
			// q = 6
			nread = neighborList[n+5*Np];
			fq = Aq[nread];
			nA += fq;
			fq = Bq[nread];
			nB += fq;

			// save the number densities
			Den[n] = nA;
			Den[Np+n] = nB;
            // save the phase field
			Phi[n] = nA-nB; 	
		}
	}
}

__global__  void dvc_ScaLBL_D3Q7_AAeven_GreyscaleColorDensity(double *Aq, double *Bq, double *Den, double *Phi, int start, int finish, int Np){
	int n;
	double fq,nA,nB;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {
			// compute number density for each component
			// q=0
			fq = Aq[n];
			nA = fq;
			fq = Bq[n];
			nB = fq;
			
			// q=1
			fq = Aq[2*Np+n];
			nA += fq;
			fq = Bq[2*Np+n];
			nB += fq;

			// q=2
			fq = Aq[1*Np+n];
			nA += fq;
			fq = Bq[1*Np+n];
			nB += fq;

			// q=3
			fq = Aq[4*Np+n];
			nA += fq;
			fq = Bq[4*Np+n];
			nB += fq;

			// q = 4
			fq = Aq[3*Np+n];
			nA += fq;
			fq = Bq[3*Np+n];
			nB += fq;
			
			// q=5
			fq = Aq[6*Np+n];
			nA += fq;
			fq = Bq[6*Np+n];
			nB += fq;
			
			// q = 6
			fq = Aq[5*Np+n];
			nA += fq;
			fq = Bq[5*Np+n];
			nB += fq;

			// save the number densities
			Den[n] = nA;
			Den[Np+n] = nB;
            // save the phase field
			Phi[n] = nA-nB; 	
		}
	}
}

__global__  void dvc_ScaLBL_D3Q7_AAodd_GreyscaleColorPhi(int *neighborList, double *Cq, double *Den, double *Phi, int start, int finish, int Np){
	int n,nread;
	double fq,phi;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {
			//..........Compute the number density for each component ............
			// q=0
			fq = Cq[n];
			phi = fq;
			
			// q=1
			nread = neighborList[n]; 
			fq = Cq[nread];
			phi += fq;
			
			// q=2
			nread = neighborList[n+Np]; 
			fq = Cq[nread];  
			phi += fq;
			
			// q=3
			nread = neighborList[n+2*Np]; 
			fq = Cq[nread];
			phi += fq;
			
			// q = 4
			nread = neighborList[n+3*Np]; 
			fq = Cq[nread];
			phi += fq;

			// q=5
			nread = neighborList[n+4*Np];
			fq = Cq[nread];
			phi += fq;
			
			// q = 6
			nread = neighborList[n+5*Np];
			fq = Cq[nread];
			phi += fq;

			// save the number densities
			Den[0*Np+n] = 0.5*(1.0+phi);
			Den[1*Np+n] = 0.5*(1.0-phi);
            // save the phase field
			Phi[n] = phi; 	
		}
	}
}

__global__  void dvc_ScaLBL_D3Q7_AAeven_GreyscaleColorPhi(double *Cq, double *Den, double *Phi, int start, int finish, int Np){
	int n;
	double fq,phi;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {
			// compute number density for each component
			// q=0
			fq = Cq[n];
			phi = fq;
			
			// q=1
			fq = Cq[2*Np+n];
			phi += fq;

			// q=2
			fq = Cq[1*Np+n];
			phi += fq;

			// q=3
			fq = Cq[4*Np+n];
			phi += fq;

			// q = 4
			fq = Cq[3*Np+n];
			phi += fq;
			
			// q=5
			fq = Cq[6*Np+n];
			phi += fq;
			
			// q = 6
			fq = Cq[5*Np+n];
			phi += fq;

			// save the number densities
			Den[0*Np+n] = 0.5*(1.0+phi);
			Den[1*Np+n] = 0.5*(1.0-phi);
            // save the phase field
			Phi[n] = phi; 	
		}
	}
}

__global__ void dvc_ScaLBL_D3Q19_GreyscaleColor_Gradient(int *neighborList, double *Den, double *DenGrad, int start, int finish, int Np){

	int n,nn;
	// distributions
	double m1,m2,m3,m4,m5,m6,m7,m8,m9;
	double m10,m11,m12,m13,m14,m15,m16,m17,m18;
	double nx,ny,nz;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {
//			nn = neighborList[n+Np]%Np;
//			m1 = Den[nn]*int(n!=nn);
//			nn = neighborList[n]%Np;
//			m2 = Den[nn]*int(n!=nn);
//			nn = neighborList[n+3*Np]%Np;
//			m3 = Den[nn]*int(n!=nn);
//			nn = neighborList[n+2*Np]%Np;
//			m4 = Den[nn]*int(n!=nn);		
//			nn = neighborList[n+5*Np]%Np;
//			m5 = Den[nn]*int(n!=nn);
//			nn = neighborList[n+4*Np]%Np;
//			m6 = Den[nn]*int(n!=nn);		
//			nn = neighborList[n+7*Np]%Np;
//			m7 = Den[nn]*int(n!=nn);
//			nn = neighborList[n+6*Np]%Np;
//			m8 = Den[nn]*int(n!=nn);		
//			nn = neighborList[n+9*Np]%Np;
//			m9 = Den[nn]*int(n!=nn);
//			nn = neighborList[n+8*Np]%Np;
//			m10 = Den[nn]*int(n!=nn);		
//			nn = neighborList[n+11*Np]%Np;
//			m11 = Den[nn]*int(n!=nn);
//			nn = neighborList[n+10*Np]%Np;
//			m12 = Den[nn]*int(n!=nn);		
//			nn = neighborList[n+13*Np]%Np;
//			m13 = Den[nn]*int(n!=nn);
//			nn = neighborList[n+12*Np]%Np;
//			m14 = Den[nn]*int(n!=nn);		
//			nn = neighborList[n+15*Np]%Np;
//			m15 = Den[nn]*int(n!=nn);
//			nn = neighborList[n+14*Np]%Np;
//			m16 = Den[nn]*int(n!=nn);		
//			nn = neighborList[n+17*Np]%Np;
//			m17 = Den[nn]*int(n!=nn);
//			nn = neighborList[n+16*Np]%Np;
//			m18 = Den[nn]*int(n!=nn);					

			nn = neighborList[n+Np]%Np;
			m1 = Den[nn];
			nn = neighborList[n]%Np;
			m2 = Den[nn];
			nn = neighborList[n+3*Np]%Np;
			m3 = Den[nn];
			nn = neighborList[n+2*Np]%Np;
			m4 = Den[nn];		
			nn = neighborList[n+5*Np]%Np;
			m5 = Den[nn];
			nn = neighborList[n+4*Np]%Np;
			m6 = Den[nn];		
			nn = neighborList[n+7*Np]%Np;
			m7 = Den[nn];
			nn = neighborList[n+6*Np]%Np;
			m8 = Den[nn];		
			nn = neighborList[n+9*Np]%Np;
			m9 = Den[nn];
			nn = neighborList[n+8*Np]%Np;
			m10 = Den[nn];		
			nn = neighborList[n+11*Np]%Np;
			m11 = Den[nn];
			nn = neighborList[n+10*Np]%Np;
			m12 = Den[nn];		
			nn = neighborList[n+13*Np]%Np;
			m13 = Den[nn];
			nn = neighborList[n+12*Np]%Np;
			m14 = Den[nn];		
			nn = neighborList[n+15*Np]%Np;
			m15 = Den[nn];
			nn = neighborList[n+14*Np]%Np;
			m16 = Den[nn];		
			nn = neighborList[n+17*Np]%Np;
			m17 = Den[nn];
			nn = neighborList[n+16*Np]%Np;
			m18 = Den[nn];					

			//............Compute the Color Gradient...................................
			nx = 1.f/6.f*(m1-m2+0.5*(m7-m8+m9-m10+m11-m12+m13-m14));
			ny = 1.f/6.f*(m3-m4+0.5*(m7-m8-m9+m10+m15-m16+m17-m18));
			nz = 1.f/6.f*(m5-m6+0.5*(m11-m12-m13+m14+m15-m16-m17+m18));
			
			DenGrad[n] = nx;
			DenGrad[Np+n] = ny;
			DenGrad[2*Np+n] = nz;
		}
	}
}

__global__ void dvc_ScaLBL_D3Q19_GreyscaleColor_Laplacian(int *neighborList, double *Den, double *DenLap, int start, int finish, int Np){

	int n,nn;
	// distributions
	double m1,m2,m3,m4,m5,m6,m7,m8,m9;
	double m10,m11,m12,m13,m14,m15,m16,m17,m18;
	double lap;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {
//			nn = neighborList[n+Np]%Np;
//			m1 = Den[nn]*int(n!=nn);
//			nn = neighborList[n]%Np;
//			m2 = Den[nn]*int(n!=nn);
//			nn = neighborList[n+3*Np]%Np;
//			m3 = Den[nn]*int(n!=nn);
//			nn = neighborList[n+2*Np]%Np;
//			m4 = Den[nn]*int(n!=nn);		
//			nn = neighborList[n+5*Np]%Np;
//			m5 = Den[nn]*int(n!=nn);
//			nn = neighborList[n+4*Np]%Np;
//			m6 = Den[nn]*int(n!=nn);		
//			nn = neighborList[n+7*Np]%Np;
//			m7 = Den[nn]*int(n!=nn);
//			nn = neighborList[n+6*Np]%Np;
//			m8 = Den[nn]*int(n!=nn);		
//			nn = neighborList[n+9*Np]%Np;
//			m9 = Den[nn]*int(n!=nn);
//			nn = neighborList[n+8*Np]%Np;
//			m10 = Den[nn]*int(n!=nn);		
//			nn = neighborList[n+11*Np]%Np;
//			m11 = Den[nn]*int(n!=nn);
//			nn = neighborList[n+10*Np]%Np;
//			m12 = Den[nn]*int(n!=nn);		
//			nn = neighborList[n+13*Np]%Np;
//			m13 = Den[nn]*int(n!=nn);
//			nn = neighborList[n+12*Np]%Np;
//			m14 = Den[nn]*int(n!=nn);		
//			nn = neighborList[n+15*Np]%Np;
//			m15 = Den[nn]*int(n!=nn);
//			nn = neighborList[n+14*Np]%Np;
//			m16 = Den[nn]*int(n!=nn);		
//			nn = neighborList[n+17*Np]%Np;
//			m17 = Den[nn]*int(n!=nn);
//			nn = neighborList[n+16*Np]%Np;
//			m18 = Den[nn]*int(n!=nn);					
			
			nn = neighborList[n+Np]%Np;
			m1 = Den[nn];
			nn = neighborList[n]%Np;
			m2 = Den[nn];
			nn = neighborList[n+3*Np]%Np;
			m3 = Den[nn];
			nn = neighborList[n+2*Np]%Np;
			m4 = Den[nn];		
			nn = neighborList[n+5*Np]%Np;
			m5 = Den[nn];
			nn = neighborList[n+4*Np]%Np;
			m6 = Den[nn];		
			nn = neighborList[n+7*Np]%Np;
			m7 = Den[nn];
			nn = neighborList[n+6*Np]%Np;
			m8 = Den[nn];		
			nn = neighborList[n+9*Np]%Np;
			m9 = Den[nn];
			nn = neighborList[n+8*Np]%Np;
			m10 = Den[nn];		
			nn = neighborList[n+11*Np]%Np;
			m11 = Den[nn];
			nn = neighborList[n+10*Np]%Np;
			m12 = Den[nn];		
			nn = neighborList[n+13*Np]%Np;
			m13 = Den[nn];
			nn = neighborList[n+12*Np]%Np;
			m14 = Den[nn];		
			nn = neighborList[n+15*Np]%Np;
			m15 = Den[nn];
			nn = neighborList[n+14*Np]%Np;
			m16 = Den[nn];		
			nn = neighborList[n+17*Np]%Np;
			m17 = Den[nn];
			nn = neighborList[n+16*Np]%Np;
			m18 = Den[nn];					


            lap = 1.f/3.f*(m1+m2+m3+m4+m5+m6-6*Den[n]+0.5*(m7+m8+m9+m10+m11+m12+m13+m14+m15+m16+m17+m18-12*Den[n]));
			DenLap[n] = lap;
		}
	}
}

__global__  void dvc_ScaLBL_D3Q19_GreyscaleColor_PressureTensor(int *neighborList, double *Phi, double *PressTensor, double *PhiLap,
      		     double kappaA,double kappaB,double lambdaA,double lambdaB, int start, int finish, int Np){
	//**GreyscaleColor model related parameters:
	//kappaA, kappaB: characterize interfacial tension
	//lambdaA, lambdaB: characterize bulk free energy 
	//nA: concentration of liquid 1; 
	//nB: concentration of liquid 2;
	//nA = 0.5*(1+phi/chi)
	//nB = 0.5*(1-phi/chi)
	//nA+nB=1
	//chi: a scaling factor, is set to 1.0 for now.

	int nn,n;
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
	double m3,m5,m7;
    double nx,ny,nz;//Color gradient
    double nA,nB;//ELBM parameters: concentration of liquid 1 and 2
    double phi;//phase field
    double pb;//thermodynamic bulk fluid pressure
    double Lphi;//Laplacian of phase field
    double C;//squared magnitude of the gradient of phase field
    double chi = 1.0;//legacy ELBM parameter, scale the phase field; may be useful in the future;
    double kappa = 0.25*(kappaA+kappaB)/(chi*chi);//the effective surface tension coefficient
    double Pxx,Pyy,Pzz,Pxy,Pyz,Pxz;//Pressure tensor

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {

//			nn = neighborList[n+Np]%Np;
//			m1 = Phi[nn]*int(n!=nn);
//			nn = neighborList[n]%Np;
//			m2 = Phi[nn]*int(n!=nn);
//			nn = neighborList[n+3*Np]%Np;
//			m3 = Phi[nn]*int(n!=nn);
//			nn = neighborList[n+2*Np]%Np;
//			m4 = Phi[nn]*int(n!=nn);		
//			nn = neighborList[n+5*Np]%Np;
//			m5 = Phi[nn]*int(n!=nn);
//			nn = neighborList[n+4*Np]%Np;
//			m6 = Phi[nn]*int(n!=nn);		
//			nn = neighborList[n+7*Np]%Np;
//			m7 = Phi[nn]*int(n!=nn);
//			nn = neighborList[n+6*Np]%Np;
//			m8 = Phi[nn]*int(n!=nn);		
//			nn = neighborList[n+9*Np]%Np;
//			m9 = Phi[nn]*int(n!=nn);
//			nn = neighborList[n+8*Np]%Np;
//			m10 = Phi[nn]*int(n!=nn);		
//			nn = neighborList[n+11*Np]%Np;
//			m11 = Phi[nn]*int(n!=nn);
//			nn = neighborList[n+10*Np]%Np;
//			m12 = Phi[nn]*int(n!=nn);		
//			nn = neighborList[n+13*Np]%Np;
//			m13 = Phi[nn]*int(n!=nn);
//			nn = neighborList[n+12*Np]%Np;
//			m14 = Phi[nn]*int(n!=nn);		
//			nn = neighborList[n+15*Np]%Np;
//			m15 = Phi[nn]*int(n!=nn);
//			nn = neighborList[n+14*Np]%Np;
//			m16 = Phi[nn]*int(n!=nn);		
//			nn = neighborList[n+17*Np]%Np;
//			m17 = Phi[nn]*int(n!=nn);
//			nn = neighborList[n+16*Np]%Np;
//			m18 = Phi[nn]*int(n!=nn);					

			nn = neighborList[n+Np]%Np;
			m1 = Phi[nn];
			nn = neighborList[n]%Np;
			m2 = Phi[nn];
			nn = neighborList[n+3*Np]%Np;
			m3 = Phi[nn];
			nn = neighborList[n+2*Np]%Np;
			m4 = Phi[nn];		
			nn = neighborList[n+5*Np]%Np;
			m5 = Phi[nn];
			nn = neighborList[n+4*Np]%Np;
			m6 = Phi[nn];		
			nn = neighborList[n+7*Np]%Np;
			m7 = Phi[nn];
			nn = neighborList[n+6*Np]%Np;
			m8 = Phi[nn];		
			nn = neighborList[n+9*Np]%Np;
			m9 = Phi[nn];
			nn = neighborList[n+8*Np]%Np;
			m10 = Phi[nn];		
			nn = neighborList[n+11*Np]%Np;
			m11 = Phi[nn];
			nn = neighborList[n+10*Np]%Np;
			m12 = Phi[nn];		
			nn = neighborList[n+13*Np]%Np;
			m13 = Phi[nn];
			nn = neighborList[n+12*Np]%Np;
			m14 = Phi[nn];		
			nn = neighborList[n+15*Np]%Np;
			m15 = Phi[nn];
			nn = neighborList[n+14*Np]%Np;
			m16 = Phi[nn];		
			nn = neighborList[n+17*Np]%Np;
			m17 = Phi[nn];
			nn = neighborList[n+16*Np]%Np;
			m18 = Phi[nn];					

			//............Compute the Color Gradient...................................
			nx = 1.f/6.f*(m1-m2+0.5*(m7-m8+m9-m10+m11-m12+m13-m14));
			ny = 1.f/6.f*(m3-m4+0.5*(m7-m8-m9+m10+m15-m16+m17-m18));
			nz = 1.f/6.f*(m5-m6+0.5*(m11-m12-m13+m14+m15-m16-m17+m18));
			C = nx*nx+ny*ny+nz*nz;
			// Laplacian of phase field
			//Lphi = 0.3333333333333333*(m1+m2+m3+m4+m5+m6)+
			//		0.16666666666666666*(m7+m8+m9+m10+m11+m12+m13+m14+m15+m16+m17+m18) - 4.0*phi;
            phi = Phi[n];
            Lphi = 1.f/3.f*(m1+m2+m3+m4+m5+m6-6*phi+0.5*(m7+m8+m9+m10+m11+m12+m13+m14+m15+m16+m17+m18-12*phi));

			//bulk pressure p_b
			nA = 0.5*(1.0+phi/chi);
			nB = 0.5*(1.0-phi/chi);
            pb = -((1.0-nA)*(1.0-nA)*nA*nA*lambdaA)*0.5 - ((1.0-nB)*(1.0-nB)*nB*nB*lambdaB)*0.5 + 
                (nA - nB)*chi*(((0.5*nA-1.5*nA*nA+nA*nA*nA)*lambdaA)/chi - ((0.5*nB-1.5*nB*nB+nB*nB*nB)*lambdaB)/chi);

			//Pressure tensors
			Pxx=pb-kappa*phi*Lphi-0.5*kappa*C + kappa*nx*nx ;
			Pyy=pb-kappa*phi*Lphi-0.5*kappa*C + kappa*ny*ny ;
			Pzz=pb-kappa*phi*Lphi-0.5*kappa*C + kappa*nz*nz ;
			Pxy= kappa*nx*ny;
			Pyz= kappa*ny*nz;
			Pxz= kappa*nx*nz;

			//...Store the Pressure Tensors....................
			PressTensor[n+0*Np] = Pxx;
			PressTensor[n+1*Np] = Pyy;
			PressTensor[n+2*Np] = Pzz;
			PressTensor[n+3*Np] = Pxy;
			PressTensor[n+4*Np] = Pyz;
			PressTensor[n+5*Np] = Pxz;
			//...............................................

			//...Store the Laplacian of phase field....................
			PhiLap[n]=Lphi;
		}
	}
}


extern "C" void ScaLBL_D3Q19_AAeven_GreyscaleColor(double *dist, double *Aq, double *Bq, double *Den,
                double *DenGradA, double *DenGradB, double *SolidForce, int start, int finish, int Np,
                double tauA,double tauB,double tauA_eff,double tauB_eff,double rhoA,double rhoB,double Gsc, double Gx, double Gy, double Gz,
                double *Poros,double *Perm, double *Velocity,double *Pressure){

    dvc_ScaLBL_D3Q19_AAeven_GreyscaleColor<<<NBLOCKS,NTHREADS >>>(dist, Aq, Bq, Den, DenGradA, DenGradB, SolidForce, start, finish, Np,
                                                                 tauA, tauB, tauA_eff, tauB_eff, rhoA, rhoB, Gsc, Gx, Gy, Gz, Poros, Perm, Velocity, Pressure);

    cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_GreyscaleColor: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_AAodd_GreyscaleColor(int *neighborList, double *dist, double *Aq, double *Bq, double *Den,
                double *DenGradA, double *DenGradB, double *SolidForce, int start, int finish, int Np,
                double tauA,double tauB,double tauA_eff,double tauB_eff,double rhoA,double rhoB,double Gsc, double Gx, double Gy, double Gz,
                double *Poros,double *Perm, double *Velocity,double *Pressure){

    dvc_ScaLBL_D3Q19_AAodd_GreyscaleColor<<<NBLOCKS,NTHREADS >>>(neighborList, dist, Aq, Bq, Den, DenGradA, DenGradB, SolidForce, start, finish, Np,
                                                                 tauA, tauB, tauA_eff, tauB_eff, rhoA, rhoB, Gsc, Gx, Gy, Gz, Poros, Perm, Velocity, Pressure);

    cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAodd_GreyscaleColor: %s \n",cudaGetErrorString(err));
	}
}

//extern "C" void ScaLBL_D3Q19_AAeven_GreyscaleColorChem(double *dist, double *Aq, double *Bq, double *Den,double *SolidForce, int start, int finish, int Np,
//                double tauA,double tauB,double tauA_eff,double tauB_eff,double rhoA,double rhoB,double gamma,double kappaA,double kappaB,double lambdaA,double lambdaB,
//                double Gx, double Gy, double Gz,
//                double *Poros,double *Perm, double *Velocity,double *Pressure,double *PressureGrad,double *PressTensorGrad,double *PhiLap){
//
//    dvc_ScaLBL_D3Q19_AAeven_GreyscaleColorChem<<<NBLOCKS,NTHREADS >>>(dist, Aq, Bq, Den, SolidForce, start, finish, Np,
//                                                                 tauA, tauB, tauA_eff, tauB_eff, rhoA, rhoB, gamma,kappaA,kappaB,lambdaA,lambdaB, Gx, Gy, Gz, Poros, Perm, Velocity, Pressure,PressureGrad,PressTensorGrad,PhiLap);
//
//    cudaError_t err = cudaGetLastError();
//	if (cudaSuccess != err){
//		printf("CUDA error in ScaLBL_D3Q19_AAeven_GreyscaleColorChem: %s \n",cudaGetErrorString(err));
//	}
//}
//
//extern "C" void ScaLBL_D3Q19_AAodd_GreyscaleColorChem(int *neighborList, double *dist, double *Aq, double *Bq, double *Den,double *SolidForce, int start, int finish, int Np,
//                double tauA,double tauB,double tauA_eff,double tauB_eff,double rhoA,double rhoB,double gamma,double kappaA,double kappaB,double lambdaA,double lambdaB,
//                double Gx, double Gy, double Gz,
//                double *Poros,double *Perm, double *Velocity,double *Pressure,double *PressureGrad,double *PressTensorGrad,double *PhiLap){
//
//    dvc_ScaLBL_D3Q19_AAodd_GreyscaleColorChem<<<NBLOCKS,NTHREADS >>>(neighborList, dist, Aq, Bq, Den, SolidForce, start, finish, Np,
//                                                                 tauA, tauB, tauA_eff, tauB_eff, rhoA, rhoB, gamma,kappaA,kappaB,lambdaA,lambdaB, Gx, Gy, Gz, 
//                                                                 Poros, Perm, Velocity, Pressure,PressureGrad,PressTensorGrad,PhiLap);
//
//    cudaError_t err = cudaGetLastError();
//	if (cudaSuccess != err){
//		printf("CUDA error in ScaLBL_D3Q19_AAodd_GreyscaleColorChem: %s \n",cudaGetErrorString(err));
//	}
//}

extern "C" void ScaLBL_D3Q19_AAeven_GreyscaleColorChem(double *dist, double *Cq, double *Phi, double *Den,double *SolidForce, int start, int finish, int Np,
                double tauA,double tauB,double tauA_eff,double tauB_eff,double rhoA,double rhoB,double gamma,double kappaA,double kappaB,double lambdaA,double lambdaB,
                double Gx, double Gy, double Gz,
                double *Poros,double *Perm, double *Velocity,double *Pressure,double *PressureGrad,double *PressTensorGrad,double *PhiLap){

    dvc_ScaLBL_D3Q19_AAeven_GreyscaleColorChem<<<NBLOCKS,NTHREADS >>>(dist, Cq, Phi, Den, SolidForce, start, finish, Np,
                                                                 tauA, tauB, tauA_eff, tauB_eff, rhoA, rhoB, gamma,kappaA,kappaB,lambdaA,lambdaB, Gx, Gy, Gz, Poros, Perm, Velocity, Pressure,PressureGrad,PressTensorGrad,PhiLap);

    cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_GreyscaleColorChem: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_AAodd_GreyscaleColorChem(int *neighborList, double *dist, double *Cq, double *Phi, double *Den,double *SolidForce, int start, int finish, int Np,
                double tauA,double tauB,double tauA_eff,double tauB_eff,double rhoA,double rhoB,double gamma,double kappaA,double kappaB,double lambdaA,double lambdaB,
                double Gx, double Gy, double Gz,
                double *Poros,double *Perm, double *Velocity,double *Pressure,double *PressureGrad,double *PressTensorGrad,double *PhiLap){

    dvc_ScaLBL_D3Q19_AAodd_GreyscaleColorChem<<<NBLOCKS,NTHREADS >>>(neighborList, dist, Cq, Phi, Den, SolidForce, start, finish, Np,
                                                                 tauA, tauB, tauA_eff, tauB_eff, rhoA, rhoB, gamma,kappaA,kappaB,lambdaA,lambdaB, Gx, Gy, Gz, 
                                                                 Poros, Perm, Velocity, Pressure,PressureGrad,PressTensorGrad,PhiLap);

    cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAodd_GreyscaleColorChem: %s \n",cudaGetErrorString(err));
	}
}

//extern "C" void ScaLBL_D3Q7_GreyColorIMRT_Init(double *Den, double *Aq, double *Bq, double *Phi, int start, int finish, int Np){
//	dvc_ScaLBL_D3Q7_GreyColorIMRT_Init<<<NBLOCKS,NTHREADS >>>(Den, Aq, Bq, Phi, start, finish, Np);
//	cudaError_t err = cudaGetLastError();
//	if (cudaSuccess != err){
//		printf("CUDA error in ScaLBL_D3Q7_GreyColorIMRT_Init: %s \n",cudaGetErrorString(err));
//	}
//}

extern "C" void ScaLBL_D3Q7_GreyColorIMRT_Init(double *Den, double *Cq, double *PhiLap, double gamma, double kappaA, double kappaB, double lambdaA, double lambdaB, int start, int finish, int Np){
	dvc_ScaLBL_D3Q7_GreyColorIMRT_Init<<<NBLOCKS,NTHREADS >>>(Den, Cq, PhiLap,gamma,kappaA,kappaB,lambdaA,lambdaB, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_GreyColorIMRT_Init: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_GreyColorIMRT_Init(double *dist, double *Den, double rhoA, double rhoB, int Np){
	dvc_ScaLBL_D3Q19_GreyColorIMRT_Init<<<NBLOCKS,NTHREADS >>>(dist,Den,rhoA,rhoB,Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_GreyColorIMRT_Init: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q7_AAodd_GreyscaleColorDensity(int *NeighborList, double *Aq, double *Bq, double *Den, double *Phi, int start, int finish, int Np){

	dvc_ScaLBL_D3Q7_AAodd_GreyscaleColorDensity<<<NBLOCKS,NTHREADS >>>(NeighborList, Aq, Bq, Den, Phi, start, finish, Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_AAodd_GreyscaleColorDensity: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q7_AAeven_GreyscaleColorDensity(double *Aq, double *Bq, double *Den, double *Phi, int start, int finish, int Np){

	dvc_ScaLBL_D3Q7_AAeven_GreyscaleColorDensity<<<NBLOCKS,NTHREADS >>>(Aq, Bq, Den, Phi, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_AAeven_GreyscaleColorDensity: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q7_AAodd_GreyscaleColorPhi(int *NeighborList, double *Cq, double *Den, double *Phi, int start, int finish, int Np){

	dvc_ScaLBL_D3Q7_AAodd_GreyscaleColorPhi<<<NBLOCKS,NTHREADS >>>(NeighborList, Cq, Den, Phi, start, finish, Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_AAodd_GreyscaleColorPhi: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q7_AAeven_GreyscaleColorPhi(double *Cq, double *Den, double *Phi, int start, int finish, int Np){

	dvc_ScaLBL_D3Q7_AAeven_GreyscaleColorPhi<<<NBLOCKS,NTHREADS >>>(Cq, Den, Phi, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_AAeven_GreyscaleColorPhi: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_GreyscaleColor_Gradient(int *neighborList, double *Den, double *DenGrad, int start, int finish, int Np){

	dvc_ScaLBL_D3Q19_GreyscaleColor_Gradient<<<NBLOCKS,NTHREADS >>>(neighborList, Den, DenGrad, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_GreyscaleColor_Gradient: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_GreyscaleColor_Laplacian(int *neighborList, double *Den, double *DenLap, int start, int finish, int Np){
	dvc_ScaLBL_D3Q19_GreyscaleColor_Laplacian<<<NBLOCKS,NTHREADS >>>(neighborList, Den, DenLap, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_GreyscaleColor_Laplacian: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_GreyscaleColor_Pressure(double *dist, double *Den, double *Porosity,double *Velocity,
                double *Pressure, double rhoA,double rhoB, int Np){

	dvc_ScaLBL_D3Q19_GreyscaleColor_Pressure<<<NBLOCKS,NTHREADS >>>(dist, Den, Porosity, Velocity, Pressure, rhoA, rhoB, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_GreyscaleColor_Pressure: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_GreyscaleColor_PressureTensor(int *neighborList, double *Phi, double *PressTensor, double *PhiLap,
      		     double kappaA,double kappaB,double lambdaA,double lambdaB, int start, int finish, int Np){
	dvc_ScaLBL_D3Q19_GreyscaleColor_PressureTensor<<<NBLOCKS,NTHREADS >>>(neighborList,Phi,PressTensor,PhiLap,kappaA,kappaB,lambdaA,lambdaB,start,finish,Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_GreyscaleColor_PressureTensor: %s \n",cudaGetErrorString(err));
	}
}
