#include <stdio.h>

#define NBLOCKS 1024
#define NTHREADS 256

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
            Gff_x = -Gsc*nA*nB_gradx*int(phi>0.0)-Gsc*nB*nA_gradx*int(phi<0.0);
            Gff_y = -Gsc*nA*nB_grady*int(phi>0.0)-Gsc*nB*nA_grady*int(phi<0.0);
            Gff_z = -Gsc*nA*nB_gradz*int(phi>0.0)-Gsc*nB*nA_gradz*int(phi<0.0);
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
            Gff_x = -Gsc*nA*nB_gradx*int(phi>0.0)-Gsc*nB*nA_gradx*int(phi<0.0);
            Gff_y = -Gsc*nA*nB_grady*int(phi>0.0)-Gsc*nB*nA_grady*int(phi<0.0);
            Gff_z = -Gsc*nA*nB_gradz*int(phi>0.0)-Gsc*nB*nA_gradz*int(phi<0.0);
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

__global__ void dvc_ScaLBL_D3Q7_GreyColorIMRT_Init(double *Den, double *Aq, double *Bq, int start, int finish, int Np){
	int idx;
    double nA,nB;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		idx =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (idx<finish) {
            nA = Den[idx];
            nB = Den[idx+Np];

			Aq[idx]=0.3333333333333333*nA;
			Aq[Np+idx]=0.1111111111111111*nA;
			Aq[2*Np+idx]=0.1111111111111111*nA;
			Aq[3*Np+idx]=0.1111111111111111*nA;
			Aq[4*Np+idx]=0.1111111111111111*nA;
			Aq[5*Np+idx]=0.1111111111111111*nA;
			Aq[6*Np+idx]=0.1111111111111111*nA;

			Bq[idx]=0.3333333333333333*nB;
			Bq[Np+idx]=0.1111111111111111*nB;
			Bq[2*Np+idx]=0.1111111111111111*nB;
			Bq[3*Np+idx]=0.1111111111111111*nB;
			Bq[4*Np+idx]=0.1111111111111111*nB;
			Bq[5*Np+idx]=0.1111111111111111*nB;
			Bq[6*Np+idx]=0.1111111111111111*nB;
		}
	}
}

__global__  void dvc_ScaLBL_D3Q7_AAodd_GreyscaleColorDensity(int *neighborList, double *Aq, double *Bq, double *Den, int start, int finish, int Np){
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
		}
	}
}

__global__  void dvc_ScaLBL_D3Q7_AAeven_GreyscaleColorDensity(double *Aq, double *Bq, double *Den, int start, int finish, int Np){
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
			nn = neighborList[n+Np]%Np;
			m1 = Den[nn]*int(n!=nn);
			nn = neighborList[n]%Np;
			m2 = Den[nn]*int(n!=nn);
			nn = neighborList[n+3*Np]%Np;
			m3 = Den[nn]*int(n!=nn);
			nn = neighborList[n+2*Np]%Np;
			m4 = Den[nn]*int(n!=nn);		
			nn = neighborList[n+5*Np]%Np;
			m5 = Den[nn]*int(n!=nn);
			nn = neighborList[n+4*Np]%Np;
			m6 = Den[nn]*int(n!=nn);		
			nn = neighborList[n+7*Np]%Np;
			m7 = Den[nn]*int(n!=nn);
			nn = neighborList[n+6*Np]%Np;
			m8 = Den[nn]*int(n!=nn);		
			nn = neighborList[n+9*Np]%Np;
			m9 = Den[nn]*int(n!=nn);
			nn = neighborList[n+8*Np]%Np;
			m10 = Den[nn]*int(n!=nn);		
			nn = neighborList[n+11*Np]%Np;
			m11 = Den[nn]*int(n!=nn);
			nn = neighborList[n+10*Np]%Np;
			m12 = Den[nn]*int(n!=nn);		
			nn = neighborList[n+13*Np]%Np;
			m13 = Den[nn]*int(n!=nn);
			nn = neighborList[n+12*Np]%Np;
			m14 = Den[nn]*int(n!=nn);		
			nn = neighborList[n+15*Np]%Np;
			m15 = Den[nn]*int(n!=nn);
			nn = neighborList[n+14*Np]%Np;
			m16 = Den[nn]*int(n!=nn);		
			nn = neighborList[n+17*Np]%Np;
			m17 = Den[nn]*int(n!=nn);
			nn = neighborList[n+16*Np]%Np;
			m18 = Den[nn]*int(n!=nn);					
			
			//............Compute the Color Gradient...................................
			nx = 1.f/18.f*(m1-m2+0.5*(m7-m8+m9-m10+m11-m12+m13-m14));
			ny = 1.f/18.f*(m3-m4+0.5*(m7-m8-m9+m10+m15-m16+m17-m18));
			nz = 1.f/18.f*(m5-m6+0.5*(m11-m12-m13+m14+m15-m16-m17+m18));
			
			DenGrad[n] = nx;
			DenGrad[Np+n] = ny;
			DenGrad[2*Np+n] = nz;
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

extern "C" void ScaLBL_D3Q7_GreyColorIMRT_Init(double *Den, double *Aq, double *Bq, int start, int finish, int Np){
	dvc_ScaLBL_D3Q7_GreyColorIMRT_Init<<<NBLOCKS,NTHREADS >>>(Den, Aq, Bq, start, finish, Np);
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

extern "C" void ScaLBL_D3Q7_AAodd_GreyscaleColorDensity(int *NeighborList, double *Aq, double *Bq, double *Den, int start, int finish, int Np){

	dvc_ScaLBL_D3Q7_AAodd_GreyscaleColorDensity<<<NBLOCKS,NTHREADS >>>(NeighborList, Aq, Bq, Den, start, finish, Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_AAodd_GreyscaleColorDensity: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q7_AAeven_GreyscaleColorDensity(double *Aq, double *Bq, double *Den, int start, int finish, int Np){

	dvc_ScaLBL_D3Q7_AAeven_GreyscaleColorDensity<<<NBLOCKS,NTHREADS >>>(Aq, Bq, Den, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_AAeven_GreyscaleColorDensity: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_GreyscaleColor_Gradient(int *neighborList, double *Den, double *DenGrad, int start, int finish, int Np){

	dvc_ScaLBL_D3Q19_GreyscaleColor_Gradient<<<NBLOCKS,NTHREADS >>>(neighborList, Den, DenGrad, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_GreyscaleColor_Gradient: %s \n",cudaGetErrorString(err));
	}
}

