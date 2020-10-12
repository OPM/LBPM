/*
  Copyright 2020 Equinor ASA                                                          
                                                                                      
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

extern "C" void ScaLBL_D3Q19_AAeven_Greyscale(double *dist, int start, int finish, int Np, double rlx, double rlx_eff, double Gx, double Gy, double Gz,
                                              double *Poros,double *Perm, double *Velocity, double *Pressure){
	int n;
	// conserved momemnts
	double rho,vx,vy,vz,v_mag;
    double ux,uy,uz,u_mag;
    double pressure;
    //double uu;
	// non-conserved moments
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18;
    double GeoFun;//geometric function from Guo's PRE 66, 036304 (2002)
    double porosity;
    double perm;//voxel permeability
    double c0, c1; //Guo's model parameters
    double mu_eff = (1.0/rlx_eff-0.5)/3.0;//kinematic viscosity
    double Fx, Fy, Fz;//The total body force including Brinkman force and user-specified (Gx,Gy,Gz)

	for (int n=start; n<finish; n++){
		// q=0
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f5 = dist[6*Np+n];
		f6 = dist[5*Np+n];
		f7 = dist[8*Np+n];
		f8 = dist[7*Np+n];
		f9 = dist[10*Np+n];
		f10 = dist[9*Np+n];
		f11 = dist[12*Np+n];
		f12 = dist[11*Np+n];
		f13 = dist[14*Np+n];
		f14 = dist[13*Np+n];
		f15 = dist[16*Np+n];
		f16 = dist[15*Np+n];
		f17 = dist[18*Np+n];
		f18 = dist[17*Np+n];
        
        porosity = Poros[n];
        perm = Perm[n];

        c0 = 0.5*(1.0+porosity*0.5*mu_eff/perm);
        if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
        GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
        c1 = porosity*0.5*GeoFun/sqrt(perm);
        if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes

		rho = f0+f2+f1+f4+f3+f6+f5+f8+f7+f10+f9+f12+f11+f14+f13+f16+f15+f18+f17;
        pressure = rho/porosity/3.0;
		vx = (f1-f2+f7-f8+f9-f10+f11-f12+f13-f14)/rho+0.5*porosity*Gx;
		vy = (f3-f4+f7-f8-f9+f10+f15-f16+f17-f18)/rho+0.5*porosity*Gy;
		vz = (f5-f6+f11-f12-f13+f14+f15-f16-f17+f18)/rho+0.5*porosity*Gz;
        v_mag=sqrt(vx*vx+vy*vy+vz*vz);
        ux = vx/(c0+sqrt(c0*c0+c1*v_mag));
        uy = vy/(c0+sqrt(c0*c0+c1*v_mag));
        uz = vz/(c0+sqrt(c0*c0+c1*v_mag));
        u_mag=sqrt(ux*ux+uy*uy+uz*uz);

        //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
        Fx = -porosity*mu_eff/perm*ux - porosity*GeoFun/sqrt(perm)*u_mag*ux + porosity*Gx;
        Fy = -porosity*mu_eff/perm*uy - porosity*GeoFun/sqrt(perm)*u_mag*uy + porosity*Gy;
        Fz = -porosity*mu_eff/perm*uz - porosity*GeoFun/sqrt(perm)*u_mag*uz + porosity*Gz;
        if (porosity==1.0){
            Fx=Gx;
            Fy=Gy;
            Fz=Gz;
        }

        //------------------------ BGK collison where body force has higher-order terms ----------------------------------------------------------//
//		// q=0
//		dist[n] = f0*(1.0-rlx)+ rlx*0.3333333333333333*rho*(1. - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//                  + 0.3333333333333333*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(0. - (3.*uy)/porosity) + Fz*(0. - (3.*uz)/porosity));
//
//		// q = 1
//		dist[1*Np+n] = f1*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 + 3.*ux + (4.5*ux*ux)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//            +0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(3. + (6.*ux)/porosity) + Fy*(0. - (3.*uy)/porosity) + Fz*(0. - (3.*uz)/porosity));
//
//		// q=2
//		dist[2*Np+n] = f2*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 - 3.*ux + (4.5*ux*ux)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//            +0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(-3. + (6.*ux)/porosity) + Fy*(0. - (3.*uy)/porosity) + Fz*(0. - (3.*uz)/porosity));
//
//		// q = 3
//		dist[3*Np+n] = f3*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 + 3.*uy + (4.5*uy*uy)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(3. + (6.*uy)/porosity) + Fz*(0. - (3.*uz)/porosity));
//
//		// q = 4
//		dist[4*Np+n] = f4*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 - 3.*uy + (4.5*uy*uy)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(-3. + (6.*uy)/porosity) + Fz*(0. - (3.*uz)/porosity));
//
//		// q = 5
//		dist[5*Np+n] = f5*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 + 3.*uz + (4.5*uz*uz)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(0. - (3.*uy)/porosity) + Fz*(3. + (6.*uz)/porosity));
//
//		// q = 6
//		dist[6*Np+n] = f6*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 - 3.*uz + (4.5*uz*uz)/porosity - (1.5*(ux*ux+ uy*uy + uz*uz))/porosity)
//				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(0. - (3.*uy)/porosity) + Fz*(-3. + (6.*uz)/porosity));
//
//		// q = 7
//		dist[7*Np+n] = f7*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux + uy) + (4.5*(ux + uy)*(ux + uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(3. - (3.*ux)/porosity + (9.*(ux + uy))/porosity) + Fy*(3. - (3.*uy)/porosity + (9.*(ux + uy))/porosity) + 
//  Fz*(0. - (3.*uz)/porosity));
//
//		// q = 8
//		dist[8*Np+n] = f8*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux - uy) + (4.5*(-ux - uy)*(-ux - uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(-3. - (3.*ux)/porosity - (9.*(-ux - uy))/porosity) + Fy*(-3. - (9.*(-ux - uy))/porosity - (3.*uy)/porosity) + 
//  Fz*(0. - (3.*uz)/porosity));
//
//		// q = 9
//		dist[9*Np+n] = f9*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux - uy) + (4.5*(ux - uy)*(ux - uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(3. - (3.*ux)/porosity + (9.*(ux - uy))/porosity) + Fy*(-3. - (9.*(ux - uy))/porosity - (3.*uy)/porosity) + 
//  Fz*(0. - (3.*uz)/porosity));
//
//		// q = 10
//		dist[10*Np+n] = f10*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux + uy) + (4.5*(-ux + uy)*(-ux + uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(-3. - (3.*ux)/porosity - (9.*(-ux + uy))/porosity) + Fy*(3. - (3.*uy)/porosity + (9.*(-ux + uy))/porosity) + 
//  Fz*(0. - (3.*uz)/porosity));
//
//		// q = 11
//		dist[11*Np+n] = f11*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux + uz) + (4.5*(ux + uz)*(ux + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(0. - (3.*uy)/porosity) + Fx*(3. - (3.*ux)/porosity + (9.*(ux + uz))/porosity) + 
//  Fz*(3. - (3.*uz)/porosity + (9.*(ux + uz))/porosity));
//
//		// q = 12
//		dist[12*Np+n] = f12*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux - uz) + (4.5*(-ux - uz)*(-ux - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(0. - (3.*uy)/porosity) + Fx*(-3. - (3.*ux)/porosity - (9.*(-ux - uz))/porosity) + 
//  Fz*(-3. - (9.*(-ux - uz))/porosity - (3.*uz)/porosity));
//
//		// q = 13
//		dist[13*Np+n] = f13*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux - uz) + (4.5*(ux - uz)*(ux - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(0. - (3.*uy)/porosity) + Fx*(3. - (3.*ux)/porosity + (9.*(ux - uz))/porosity) + 
//  Fz*(-3. - (9.*(ux - uz))/porosity - (3.*uz)/porosity));
//
//		// q= 14
//		dist[14*Np+n] = f14*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux + uz) + (4.5*(-ux + uz)*(-ux + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(0. - (3.*uy)/porosity) + Fx*(-3. - (3.*ux)/porosity - (9.*(-ux + uz))/porosity) + 
//  Fz*(3. - (3.*uz)/porosity + (9.*(-ux + uz))/porosity));
//
//		// q = 15
//		dist[15*Np+n] = f15*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(uy + uz) + (4.5*(uy + uz)*(uy + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(3. - (3.*uy)/porosity + (9.*(uy + uz))/porosity) + 
//  Fz*(3. - (3.*uz)/porosity + (9.*(uy + uz))/porosity));
//
//		// q = 16
//		dist[16*Np+n] = f16*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-uy - uz) + (4.5*(-uy - uz)*(-uy - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(-3. - (3.*uy)/porosity - (9.*(-uy - uz))/porosity) + 
//  Fz*(-3. - (9.*(-uy - uz))/porosity - (3.*uz)/porosity));
//
//		// q = 17
//		dist[17*Np+n] = f17*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(uy - uz) + (4.5*(uy - uz)*(uy - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(3. - (3.*uy)/porosity + (9.*(uy - uz))/porosity) + 
//  Fz*(-3. - (9.*(uy - uz))/porosity - (3.*uz)/porosity));
//
//		// q = 18
//		dist[18*Np+n] = f18*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-uy + uz) + (4.5*(-uy + uz)*(-uy + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(-3. - (3.*uy)/porosity - (9.*(-uy + uz))/porosity) + 
//  Fz*(3. - (3.*uz)/porosity + (9.*(-uy + uz))/porosity));
        //----------------------------------------------------------------------------------------------------------------------------------------//

        //------------------------ BGK collison where body force has NO higher-order terms ----------------------------------------------------------//
		// q=0
		dist[n] = f0*(1.0-rlx)+ rlx*0.3333333333333333*rho*(1. - (1.5*(ux*ux + uy*uy + uz*uz))/porosity);

		// q = 1
		dist[1*Np+n] = f1*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 + 3.*ux + (4.5*ux*ux)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
            +0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(3.));

		// q=2
		dist[2*Np+n] = f2*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 - 3.*ux + (4.5*ux*ux)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
            +0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(-3.));

		// q = 3
		dist[3*Np+n] = f3*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 + 3.*uy + (4.5*uy*uy)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fy*(3.));

		// q = 4
		dist[4*Np+n] = f4*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 - 3.*uy + (4.5*uy*uy)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fy*(-3.));

		// q = 5
		dist[5*Np+n] = f5*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 + 3.*uz + (4.5*uz*uz)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fz*(3.));

		// q = 6
		dist[6*Np+n] = f6*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 - 3.*uz + (4.5*uz*uz)/porosity - (1.5*(ux*ux+ uy*uy + uz*uz))/porosity)
				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fz*(-3.));

		// q = 7
		dist[7*Np+n] = f7*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux + uy) + (4.5*(ux + uy)*(ux + uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(3.) + Fy*(3.));

		// q = 8
		dist[8*Np+n] = f8*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux - uy) + (4.5*(-ux - uy)*(-ux - uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(-3.) + Fy*(-3.));

		// q = 9
		dist[9*Np+n] = f9*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux - uy) + (4.5*(ux - uy)*(ux - uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(3.) + Fy*(-3.));

		// q = 10
		dist[10*Np+n] = f10*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux + uy) + (4.5*(-ux + uy)*(-ux + uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(-3.) + Fy*(3.));

		// q = 11
		dist[11*Np+n] = f11*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux + uz) + (4.5*(ux + uz)*(ux + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(3.) + Fz*(3.));

		// q = 12
		dist[12*Np+n] = f12*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux - uz) + (4.5*(-ux - uz)*(-ux - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(-3.) + Fz*(-3.));

		// q = 13
		dist[13*Np+n] = f13*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux - uz) + (4.5*(ux - uz)*(ux - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(3.) + Fz*(-3.));

		// q= 14
		dist[14*Np+n] = f14*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux + uz) + (4.5*(-ux + uz)*(-ux + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(-3.) + Fz*(3.));

		// q = 15
		dist[15*Np+n] = f15*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(uy + uz) + (4.5*(uy + uz)*(uy + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(3.) + Fz*(3.));

		// q = 16
		dist[16*Np+n] = f16*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-uy - uz) + (4.5*(-uy - uz)*(-uy - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(-3.) + Fz*(-3.));

		// q = 17
		dist[17*Np+n] = f17*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(uy - uz) + (4.5*(uy - uz)*(uy - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(3.) + Fz*(-3.));

		// q = 18
		dist[18*Np+n] = f18*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-uy + uz) + (4.5*(-uy + uz)*(-uy + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(-3.) + Fz*(3.));
        //-------------------------------------------------------------------------------------------------------------------------------------------//
        
        //Update velocity on device
		Velocity[0*Np+n] = ux;
		Velocity[1*Np+n] = uy;
		Velocity[2*Np+n] = uz;
        //Update pressure on device
        Pressure[n] = pressure;
	}
}

extern "C" void ScaLBL_D3Q19_AAodd_Greyscale(int *neighborList, double *dist, int start, int finish, int Np, double rlx,  double rlx_eff, double Gx, double Gy, double Gz, 
                                             double *Poros,double *Perm, double *Velocity,double *Pressure){
	int n;
	// conserved momemnts
	double rho,vx,vy,vz,v_mag;
    double ux,uy,uz,u_mag;
    double pressure;
    //double uu;
	// non-conserved moments
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18;
	int nr1,nr2,nr3,nr4,nr5,nr6,nr7,nr8,nr9,nr10,nr11,nr12,nr13,nr14,nr15,nr16,nr17,nr18;
    double GeoFun;//geometric function from Guo's PRE 66, 036304 (2002)
    double porosity;
    double perm;//voxel permeability
    double c0, c1; //Guo's model parameters
    double mu_eff = (1.0/rlx_eff-0.5)/3.0;//kinematic viscosity
    double Fx, Fy, Fz;//The total body force including Brinkman force and user-specified (Gx,Gy,Gz)

	int nread;
	for (int n=start; n<finish; n++){
		
		// q=0
		f0 = dist[n];
		// q=1
		nr1 = neighborList[n]; // neighbor 2 ( > 10Np => odd part of dist)
		f1 = dist[nr1]; // reading the f1 data into register fq

		nr2 = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
		f2 = dist[nr2];  // reading the f2 data into register fq

		// q=3
		nr3 = neighborList[n+2*Np]; // neighbor 4
		f3 = dist[nr3];

		// q = 4
		nr4 = neighborList[n+3*Np]; // neighbor 3
		f4 = dist[nr4];

		// q=5
		nr5 = neighborList[n+4*Np];
		f5 = dist[nr5];

		// q = 6
		nr6 = neighborList[n+5*Np];
		f6 = dist[nr6];
		
		// q=7
		nr7 = neighborList[n+6*Np];
		f7 = dist[nr7];

		// q = 8
		nr8 = neighborList[n+7*Np];
		f8 = dist[nr8];

		// q=9
		nr9 = neighborList[n+8*Np];
		f9 = dist[nr9];

		// q = 10
		nr10 = neighborList[n+9*Np];
		f10 = dist[nr10];

		// q=11
		nr11 = neighborList[n+10*Np];
		f11 = dist[nr11];

		// q=12
		nr12 = neighborList[n+11*Np];
		f12 = dist[nr12];

		// q=13
		nr13 = neighborList[n+12*Np];
		f13 = dist[nr13];

		// q=14
		nr14 = neighborList[n+13*Np];
		f14 = dist[nr14];

		// q=15
		nr15 = neighborList[n+14*Np];
		f15 = dist[nr15];

		// q=16
		nr16 = neighborList[n+15*Np];
		f16 = dist[nr16];

		// q=17
		//fq = dist[18*Np+n];
		nr17 = neighborList[n+16*Np];
		f17 = dist[nr17];

		// q=18
		nr18 = neighborList[n+17*Np];
		f18 = dist[nr18];

        porosity = Poros[n];
        perm = Perm[n];

        c0 = 0.5*(1.0+porosity*0.5*mu_eff/perm);
        if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
        GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
        c1 = porosity*0.5*GeoFun/sqrt(perm);
        if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes

		rho = f0+f2+f1+f4+f3+f6+f5+f8+f7+f10+f9+f12+f11+f14+f13+f16+f15+f18+f17;
        pressure = rho/porosity/3.0;
		vx = (f1-f2+f7-f8+f9-f10+f11-f12+f13-f14)/rho+0.5*porosity*Gx;
		vy = (f3-f4+f7-f8-f9+f10+f15-f16+f17-f18)/rho+0.5*porosity*Gy;
		vz = (f5-f6+f11-f12-f13+f14+f15-f16-f17+f18)/rho+0.5*porosity*Gz;
        v_mag=sqrt(vx*vx+vy*vy+vz*vz);
        ux = vx/(c0+sqrt(c0*c0+c1*v_mag));
        uy = vy/(c0+sqrt(c0*c0+c1*v_mag));
        uz = vz/(c0+sqrt(c0*c0+c1*v_mag));
        u_mag=sqrt(ux*ux+uy*uy+uz*uz);

        //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
        Fx = -porosity*mu_eff/perm*ux - porosity*GeoFun/sqrt(perm)*u_mag*ux + porosity*Gx;
        Fy = -porosity*mu_eff/perm*uy - porosity*GeoFun/sqrt(perm)*u_mag*uy + porosity*Gy;
        Fz = -porosity*mu_eff/perm*uz - porosity*GeoFun/sqrt(perm)*u_mag*uz + porosity*Gz;
        if (porosity==1.0){
            Fx=Gx;
            Fy=Gy;
            Fz=Gz;
        }

        //------------------------ BGK collison where body force has higher-order terms ----------------------------------------------------------//
//		// q=0
//		dist[n] = f0*(1.0-rlx) + rlx*0.3333333333333333*rho*(1. - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//                  + 0.3333333333333333*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(0. - (3.*uy)/porosity) + Fz*(0. - (3.*uz)/porosity));
//
//		// q = 1
//		dist[nr2] = f1*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 + 3.*ux + (4.5*ux*ux)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//            +0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(3. + (6.*ux)/porosity) + Fy*(0. - (3.*uy)/porosity) + Fz*(0. - (3.*uz)/porosity));
//
//		// q=2
//		dist[nr1] = f2*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 - 3.*ux + (4.5*ux*ux)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//            +0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(-3. + (6.*ux)/porosity) + Fy*(0. - (3.*uy)/porosity) + Fz*(0. - (3.*uz)/porosity));
//
//		// q = 3
//		dist[nr4] = f3*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 + 3.*uy + (4.5*uy*uy)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(3. + (6.*uy)/porosity) + Fz*(0. - (3.*uz)/porosity));
//
//		// q = 4
//		dist[nr3] = f4*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 - 3.*uy + (4.5*uy*uy)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)  
//				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(-3. + (6.*uy)/porosity) + Fz*(0. - (3.*uz)/porosity));
//
//		// q = 5
//		dist[nr6] = f5*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 + 3.*uz + (4.5*uz*uz)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(0. - (3.*uy)/porosity) + Fz*(3. + (6.*uz)/porosity));
//
//		// q = 6
//		dist[nr5] = f6*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 - 3.*uz + (4.5*uz*uz)/porosity - (1.5*(ux*ux+ uy*uy + uz*uz))/porosity) 
//				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(0. - (3.*uy)/porosity) + Fz*(-3. + (6.*uz)/porosity));
//
//		// q = 7
//		dist[nr8] = f7*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux + uy) + (4.5*(ux + uy)*(ux + uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(3. - (3.*ux)/porosity + (9.*(ux + uy))/porosity) + Fy*(3. - (3.*uy)/porosity + (9.*(ux + uy))/porosity) + 
//  Fz*(0. - (3.*uz)/porosity));
//
//		// q = 8
//		dist[nr7] = f8*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux - uy) + (4.5*(-ux - uy)*(-ux - uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(-3. - (3.*ux)/porosity - (9.*(-ux - uy))/porosity) + Fy*(-3. - (9.*(-ux - uy))/porosity - (3.*uy)/porosity) + 
//  Fz*(0. - (3.*uz)/porosity));
//
//		// q = 9
//		dist[nr10] = f9*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux - uy) + (4.5*(ux - uy)*(ux - uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(3. - (3.*ux)/porosity + (9.*(ux - uy))/porosity) + Fy*(-3. - (9.*(ux - uy))/porosity - (3.*uy)/porosity) + 
//  Fz*(0. - (3.*uz)/porosity));
//
//		// q = 10
//		dist[nr9] = f10*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux + uy) + (4.5*(-ux + uy)*(-ux + uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(-3. - (3.*ux)/porosity - (9.*(-ux + uy))/porosity) + Fy*(3. - (3.*uy)/porosity + (9.*(-ux + uy))/porosity) + 
//  Fz*(0. - (3.*uz)/porosity));
//
//		// q = 11
//		dist[nr12] = f11*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux + uz) + (4.5*(ux + uz)*(ux + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(0. - (3.*uy)/porosity) + Fx*(3. - (3.*ux)/porosity + (9.*(ux + uz))/porosity) + 
//  Fz*(3. - (3.*uz)/porosity + (9.*(ux + uz))/porosity));
//
//		// q = 12
//		dist[nr11] = f12*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux - uz) + (4.5*(-ux - uz)*(-ux - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(0. - (3.*uy)/porosity) + Fx*(-3. - (3.*ux)/porosity - (9.*(-ux - uz))/porosity) + 
//  Fz*(-3. - (9.*(-ux - uz))/porosity - (3.*uz)/porosity));
//
//		// q = 13
//		dist[nr14] = f13*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux - uz) + (4.5*(ux - uz)*(ux - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(0. - (3.*uy)/porosity) + Fx*(3. - (3.*ux)/porosity + (9.*(ux - uz))/porosity) + 
//  Fz*(-3. - (9.*(ux - uz))/porosity - (3.*uz)/porosity));
//
//		// q= 14
//		dist[nr13] = f14*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux + uz) + (4.5*(-ux + uz)*(-ux + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(0. - (3.*uy)/porosity) + Fx*(-3. - (3.*ux)/porosity - (9.*(-ux + uz))/porosity) + 
//  Fz*(3. - (3.*uz)/porosity + (9.*(-ux + uz))/porosity));
//
//		// q = 15
//		dist[nr16] = f15*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(uy + uz) + (4.5*(uy + uz)*(uy + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(3. - (3.*uy)/porosity + (9.*(uy + uz))/porosity) + 
//  Fz*(3. - (3.*uz)/porosity + (9.*(uy + uz))/porosity));
//
//		// q = 16
//		dist[nr15] = f16*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-uy - uz) + (4.5*(-uy - uz)*(-uy - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(-3. - (3.*uy)/porosity - (9.*(-uy - uz))/porosity) + 
//  Fz*(-3. - (9.*(-uy - uz))/porosity - (3.*uz)/porosity));
//
//		// q = 17
//		dist[nr18] = f17*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(uy - uz) + (4.5*(uy - uz)*(uy - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(3. - (3.*uy)/porosity + (9.*(uy - uz))/porosity) + 
//  Fz*(-3. - (9.*(uy - uz))/porosity - (3.*uz)/porosity));
//
//		// q = 18
//		dist[nr17] = f18*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-uy + uz) + (4.5*(-uy + uz)*(-uy + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
//				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(-3. - (3.*uy)/porosity - (9.*(-uy + uz))/porosity) + 
//  Fz*(3. - (3.*uz)/porosity + (9.*(-uy + uz))/porosity));
        //----------------------------------------------------------------------------------------------------------------------------------------//



        //------------------------ BGK collison where body force has NO higher-order terms ----------------------------------------------------------//
		// q=0
		dist[n] = f0*(1.0-rlx) + rlx*0.3333333333333333*rho*(1. - (1.5*(ux*ux + uy*uy + uz*uz))/porosity);

		// q = 1
		dist[nr2] = f1*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 + 3.*ux + (4.5*ux*ux)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
            +0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(3.));

		// q=2
		dist[nr1] = f2*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 - 3.*ux + (4.5*ux*ux)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
            +0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(-3.));

		// q = 3
		dist[nr4] = f3*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 + 3.*uy + (4.5*uy*uy)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fy*(3.));

		// q = 4
		dist[nr3] = f4*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 - 3.*uy + (4.5*uy*uy)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)  
				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fy*(-3.));

		// q = 5
		dist[nr6] = f5*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 + 3.*uz + (4.5*uz*uz)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fz*(3.));

		// q = 6
		dist[nr5] = f6*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 - 3.*uz + (4.5*uz*uz)/porosity - (1.5*(ux*ux+ uy*uy + uz*uz))/porosity) 
				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fz*(-3.));

		// q = 7
		dist[nr8] = f7*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux + uy) + (4.5*(ux + uy)*(ux + uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(3.) + Fy*(3.));

		// q = 8
		dist[nr7] = f8*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux - uy) + (4.5*(-ux - uy)*(-ux - uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(-3.) + Fy*(-3.));

		// q = 9
		dist[nr10] = f9*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux - uy) + (4.5*(ux - uy)*(ux - uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(3.) + Fy*(-3.));

		// q = 10
		dist[nr9] = f10*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux + uy) + (4.5*(-ux + uy)*(-ux + uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(-3.) + Fy*(3.));

		// q = 11
		dist[nr12] = f11*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux + uz) + (4.5*(ux + uz)*(ux + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(3.) + Fz*(3.));

		// q = 12
		dist[nr11] = f12*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux - uz) + (4.5*(-ux - uz)*(-ux - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(-3.) + Fz*(-3.));

		// q = 13
		dist[nr14] = f13*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux - uz) + (4.5*(ux - uz)*(ux - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(3.) + Fz*(-3.));

		// q= 14
		dist[nr13] = f14*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux + uz) + (4.5*(-ux + uz)*(-ux + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(-3.) + Fz*(3.));

		// q = 15
		dist[nr16] = f15*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(uy + uz) + (4.5*(uy + uz)*(uy + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(3.) + Fz*(3.));

		// q = 16
		dist[nr15] = f16*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-uy - uz) + (4.5*(-uy - uz)*(-uy - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(-3.) + Fz*(-3.));

		// q = 17
		dist[nr18] = f17*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(uy - uz) + (4.5*(uy - uz)*(uy - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(3.) + Fz*(-3.));

		// q = 18
		dist[nr17] = f18*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-uy + uz) + (4.5*(-uy + uz)*(-uy + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(-3.) + Fz*(3.));
        //-------------------------------------------------------------------------------------------------------------------------------------------//



        //Update velocity on device
		Velocity[0*Np+n] = ux;
		Velocity[1*Np+n] = uy;
		Velocity[2*Np+n] = uz;
        //Update pressure on device
        Pressure[n] = pressure;
	}
}


extern "C" void ScaLBL_D3Q19_AAeven_Greyscale_IMRT(double *dist, int start, int finish, int Np, double rlx,  double rlx_eff, double Gx, double Gy, double Gz,
                                              double *Poros,double *Perm, double *Velocity, double Den,double *Pressure){
	int n;
	double vx,vy,vz,v_mag;
    double ux,uy,uz,u_mag;
    double pressure;//defined for this incompressible model
	// conserved momemnts
	double jx,jy,jz;
	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
    double fq;
	//double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18;
    double GeoFun;//geometric function from Guo's PRE 66, 036304 (2002)
    double porosity;
    double perm;//voxel permeability
    double c0, c1; //Guo's model parameters
    double mu_eff = (1.0/rlx_eff-0.5)/3.0;//kinematic viscosity
    double Fx, Fy, Fz;//The total body force including Brinkman force and user-specified (Gx,Gy,Gz)
    double rlx_setA = rlx;
    double rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);

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

        porosity = Poros[n];
        perm = Perm[n];

        c0 = 0.5*(1.0+porosity*0.5*mu_eff/perm);
        if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
        GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
        c1 = porosity*0.5*GeoFun/sqrt(perm);
        if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes

		vx = jx/Den+0.5*porosity*Gx;
		vy = jy/Den+0.5*porosity*Gy;
		vz = jz/Den+0.5*porosity*Gz;
        v_mag=sqrt(vx*vx+vy*vy+vz*vz);
        ux = vx/(c0+sqrt(c0*c0+c1*v_mag));
        uy = vy/(c0+sqrt(c0*c0+c1*v_mag));
        uz = vz/(c0+sqrt(c0*c0+c1*v_mag));
        u_mag=sqrt(ux*ux+uy*uy+uz*uz);

        //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
        Fx = Den*(-porosity*mu_eff/perm*ux - porosity*GeoFun/sqrt(perm)*u_mag*ux + porosity*Gx);
        Fy = Den*(-porosity*mu_eff/perm*uy - porosity*GeoFun/sqrt(perm)*u_mag*uy + porosity*Gy);
        Fz = Den*(-porosity*mu_eff/perm*uz - porosity*GeoFun/sqrt(perm)*u_mag*uz + porosity*Gz);
        if (porosity==1.0){
            Fx=Den*Gx;
            Fy=Den*Gy;
            Fz=Den*Gz;
        }

        //Calculate pressure for Incompressible-MRT model
        pressure=0.5/porosity*(pressure-0.5*Den*u_mag*u_mag/porosity);

        //-------------------- IMRT collison where body force has higher-order terms -------------//
//		//..............carry out relaxation process...............................................
//		m1 = m1 + rlx_setA*((-30*Den+19*Den*(ux*ux+uy*uy+uz*uz)/porosity + 57*pressure*porosity) - m1) 
//                + (1-0.5*rlx_setA)*38*(Fx*ux+Fy*uy+Fz*uz)/porosity;
//		m2 = m2 + rlx_setA*((12*Den - 5.5*Den*(ux*ux+uy*uy+uz*uz)/porosity-27*pressure*porosity) - m2)
//                + (1-0.5*rlx_setA)*11*(-Fx*ux-Fy*uy-Fz*uz)/porosity;
//        jx = jx + Fx;
//		m4 = m4 + rlx_setB*((-0.6666666666666666*ux*Den) - m4)
//                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
//        jy = jy + Fy;
//		m6 = m6 + rlx_setB*((-0.6666666666666666*uy*Den) - m6)
//                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
//        jz = jz + Fz;
//		m8 = m8 + rlx_setB*((-0.6666666666666666*uz*Den) - m8)
//                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
//		m9 = m9 + rlx_setA*((Den*(2*ux*ux-uy*uy-uz*uz)/porosity) - m9)
//                + (1-0.5*rlx_setA)*(4*Fx*ux-2*Fy*uy-2*Fz*uz)/porosity;
//		m10 = m10 + rlx_setA*(-0.5*Den*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10)
//                  + (1-0.5*rlx_setA)*(-2*Fx*ux+Fy*uy+Fz*uz)/porosity;
//		m11 = m11 + rlx_setA*((Den*(uy*uy-uz*uz)/porosity) - m11)
//                  + (1-0.5*rlx_setA)*(2*Fy*uy-2*Fz*uz)/porosity;
//		m12 = m12 + rlx_setA*(-0.5*(Den*(uy*uy-uz*uz)/porosity)- m12)
//                  + (1-0.5*rlx_setA)*(-Fy*uy+Fz*uz)/porosity;
//		m13 = m13 + rlx_setA*((Den*ux*uy/porosity) - m13)
//                  + (1-0.5*rlx_setA)*(Fy*ux+Fx*uy)/porosity;
//		m14 = m14 + rlx_setA*((Den*uy*uz/porosity) - m14)
//                  + (1-0.5*rlx_setA)*(Fz*uy+Fy*uz)/porosity;
//		m15 = m15 + rlx_setA*((Den*ux*uz/porosity) - m15)
//                  + (1-0.5*rlx_setA)*(Fz*ux+Fx*uz)/porosity;
//		m16 = m16 + rlx_setB*( - m16);
//		m17 = m17 + rlx_setB*( - m17);
//		m18 = m18 + rlx_setB*( - m18);
//		//.......................................................................................................


        //-------------------- IMRT collison where body force has NO higher-order terms -------------//
		//..............carry out relaxation process...............................................
		m1 = m1 + rlx_setA*((-30*Den+19*Den*(ux*ux+uy*uy+uz*uz)/porosity + 57*pressure*porosity) - m1);
		m2 = m2 + rlx_setA*((12*Den - 5.5*Den*(ux*ux+uy*uy+uz*uz)/porosity-27*pressure*porosity) - m2);
        jx = jx + Fx;
		m4 = m4 + rlx_setB*((-0.6666666666666666*ux*Den) - m4)
                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
        jy = jy + Fy;
		m6 = m6 + rlx_setB*((-0.6666666666666666*uy*Den) - m6)
                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
        jz = jz + Fz;
		m8 = m8 + rlx_setB*((-0.6666666666666666*uz*Den) - m8)
                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
		m9 = m9 + rlx_setA*((Den*(2*ux*ux-uy*uy-uz*uz)/porosity) - m9);
		m10 = m10 + rlx_setA*(-0.5*Den*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10);
		m11 = m11 + rlx_setA*((Den*(uy*uy-uz*uz)/porosity) - m11);
		m12 = m12 + rlx_setA*(-0.5*(Den*(uy*uy-uz*uz)/porosity)- m12);
		m13 = m13 + rlx_setA*((Den*ux*uy/porosity) - m13);
		m14 = m14 + rlx_setA*((Den*uy*uz/porosity) - m14);
		m15 = m15 + rlx_setA*((Den*ux*uz/porosity) - m15);
		m16 = m16 + rlx_setB*( - m16);
		m17 = m17 + rlx_setB*( - m17);
		m18 = m18 + rlx_setB*( - m18);
		//.......................................................................................................


		//.................inverse transformation......................................................
		// q=0
		fq = mrt_V1*Den-mrt_V2*m1+mrt_V3*m2;
		dist[n] = fq;

		// q = 1
		fq = mrt_V1*Den-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10);
		dist[1*Np+n] = fq;

		// q=2
		fq = mrt_V1*Den-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10);
		dist[2*Np+n] = fq;

		// q = 3
		fq = mrt_V1*Den-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
		dist[3*Np+n] = fq;

		// q = 4
		fq = mrt_V1*Den-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
		dist[4*Np+n] = fq;

		// q = 5
		fq = mrt_V1*Den-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
		dist[5*Np+n] = fq;

		// q = 6
		fq = mrt_V1*Den-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
		dist[6*Np+n] = fq;

		// q = 7
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m16-m17);
		dist[7*Np+n] = fq;

		// q = 8
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m17-m16);
		dist[8*Np+n] = fq;

		// q = 9
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13+0.125*(m16+m17);
		dist[9*Np+n] = fq;

		// q = 10
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13-0.125*(m16+m17);
		dist[10*Np+n] = fq;

		// q = 11
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m18-m16);
		dist[11*Np+n] = fq;

		// q = 12
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m16-m18);
		dist[12*Np+n] = fq;

		// q = 13
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12-0.25*m15-0.125*(m16+m18);
		dist[13*Np+n] = fq;

		// q= 14
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12-0.25*m15+0.125*(m16+m18);
		dist[14*Np+n] = fq;

		// q = 15
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18);
		dist[15*Np+n] = fq;

		// q = 16
		fq =  mrt_V1*Den+mrt_V9*m1+mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17);
		dist[16*Np+n] = fq;

		// q = 17
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)-mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18);
		dist[17*Np+n] = fq;

		// q = 18
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)-mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18);
		dist[18*Np+n] = fq;
		//........................................................................

        //Update velocity on device
		Velocity[0*Np+n] = ux;
		Velocity[1*Np+n] = uy;
		Velocity[2*Np+n] = uz;
        //Update pressure on device
        Pressure[n] = pressure;
	}
}

extern "C" void ScaLBL_D3Q19_AAodd_Greyscale_IMRT(int *neighborList, double *dist, int start, int finish, int Np, double rlx,  double rlx_eff, double Gx, double Gy, double Gz, 
                                             double *Poros,double *Perm, double *Velocity, double Den,double *Pressure){
	int n, nread;
	double vx,vy,vz,v_mag;
    double ux,uy,uz,u_mag;
    double pressure;//defined for this incompressible model
	// conserved momemnts
	double jx,jy,jz;
	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
    double fq;
	//double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18;
    double GeoFun;//geometric function from Guo's PRE 66, 036304 (2002)
    double porosity;
    double perm;//voxel permeability
    double c0, c1; //Guo's model parameters
    double mu_eff = (1.0/rlx_eff-0.5)/3.0;//kinematic viscosity
    double Fx, Fy, Fz;//The total body force including Brinkman force and user-specified (Gx,Gy,Gz)
    double rlx_setA = rlx;
    double rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);

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

		//........................................................................
		//					READ THE DISTRIBUTIONS
		//		(read from opposite array due to previous swap operation)
		//........................................................................
		// q=0
		fq = dist[n];
		m1  = -30.0*fq;
		m2  = 12.0*fq;

		// q=1
		nread = neighborList[n]; // neighbor 2 ( > 10Np => odd part of dist)
		fq = dist[nread]; // reading the f1 data into register fq
        pressure = fq;
		m1 -= 11.0*fq;
		m2 -= 4.0*fq;
		jx = fq;
		m4 = -4.0*fq;
		m9 = 2.0*fq;
		m10 = -4.0*fq;

		// q=2
		nread = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
		fq = dist[nread];  // reading the f2 data into register fq
		pressure += fq;
		m1 -= 11.0*(fq);
		m2 -= 4.0*(fq);
		jx -= fq;
		m4 += 4.0*(fq);
		m9 += 2.0*(fq);
		m10 -= 4.0*(fq);

		// q=3
		nread = neighborList[n+2*Np]; // neighbor 4
		fq = dist[nread];
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
		nread = neighborList[n+3*Np]; // neighbor 3
		fq = dist[nread];
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
		nread = neighborList[n+4*Np];
		fq = dist[nread];
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
		nread = neighborList[n+5*Np];
		fq = dist[nread];
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

        porosity = Poros[n];
        perm = Perm[n];

        c0 = 0.5*(1.0+porosity*0.5*mu_eff/perm);
        if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
        GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
        c1 = porosity*0.5*GeoFun/sqrt(perm);
        if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes

		vx = jx/Den+0.5*porosity*Gx;
		vy = jy/Den+0.5*porosity*Gy;
		vz = jz/Den+0.5*porosity*Gz;
        v_mag=sqrt(vx*vx+vy*vy+vz*vz);
        ux = vx/(c0+sqrt(c0*c0+c1*v_mag));
        uy = vy/(c0+sqrt(c0*c0+c1*v_mag));
        uz = vz/(c0+sqrt(c0*c0+c1*v_mag));
        u_mag=sqrt(ux*ux+uy*uy+uz*uz);

        //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
        Fx = Den*(-porosity*mu_eff/perm*ux - porosity*GeoFun/sqrt(perm)*u_mag*ux + porosity*Gx);
        Fy = Den*(-porosity*mu_eff/perm*uy - porosity*GeoFun/sqrt(perm)*u_mag*uy + porosity*Gy);
        Fz = Den*(-porosity*mu_eff/perm*uz - porosity*GeoFun/sqrt(perm)*u_mag*uz + porosity*Gz);
        if (porosity==1.0){
            Fx=Den*Gx;
            Fy=Den*Gy;
            Fz=Den*Gz;
        }

        //Calculate pressure for Incompressible-MRT model
        pressure=0.5/porosity*(pressure-0.5*Den*u_mag*u_mag/porosity);

        //-------------------- IMRT collison where body force has higher-order terms -------------//
//		//..............carry out relaxation process...............................................
//		m1 = m1 + rlx_setA*((-30*Den+19*Den*(ux*ux+uy*uy+uz*uz)/porosity + 57*pressure*porosity) - m1) 
//                + (1-0.5*rlx_setA)*38*(Fx*ux+Fy*uy+Fz*uz)/porosity;
//		m2 = m2 + rlx_setA*((12*Den - 5.5*Den*(ux*ux+uy*uy+uz*uz)/porosity-27*pressure*porosity) - m2)
//                + (1-0.5*rlx_setA)*11*(-Fx*ux-Fy*uy-Fz*uz)/porosity;
//        jx = jx + Fx;
//		m4 = m4 + rlx_setB*((-0.6666666666666666*ux*Den) - m4)
//                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
//        jy = jy + Fy;
//		m6 = m6 + rlx_setB*((-0.6666666666666666*uy*Den) - m6)
//                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
//        jz = jz + Fz;
//		m8 = m8 + rlx_setB*((-0.6666666666666666*uz*Den) - m8)
//                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
//		m9 = m9 + rlx_setA*((Den*(2*ux*ux-uy*uy-uz*uz)/porosity) - m9)
//                + (1-0.5*rlx_setA)*(4*Fx*ux-2*Fy*uy-2*Fz*uz)/porosity;
//		m10 = m10 + rlx_setA*(-0.5*Den*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10)
//                  + (1-0.5*rlx_setA)*(-2*Fx*ux+Fy*uy+Fz*uz)/porosity;
//		m11 = m11 + rlx_setA*((Den*(uy*uy-uz*uz)/porosity) - m11)
//                  + (1-0.5*rlx_setA)*(2*Fy*uy-2*Fz*uz)/porosity;
//		m12 = m12 + rlx_setA*(-0.5*(Den*(uy*uy-uz*uz)/porosity)- m12)
//                  + (1-0.5*rlx_setA)*(-Fy*uy+Fz*uz)/porosity;
//		m13 = m13 + rlx_setA*((Den*ux*uy/porosity) - m13)
//                  + (1-0.5*rlx_setA)*(Fy*ux+Fx*uy)/porosity;
//		m14 = m14 + rlx_setA*((Den*uy*uz/porosity) - m14)
//                  + (1-0.5*rlx_setA)*(Fz*uy+Fy*uz)/porosity;
//		m15 = m15 + rlx_setA*((Den*ux*uz/porosity) - m15)
//                  + (1-0.5*rlx_setA)*(Fz*ux+Fx*uz)/porosity;
//		m16 = m16 + rlx_setB*( - m16);
//		m17 = m17 + rlx_setB*( - m17);
//		m18 = m18 + rlx_setB*( - m18);
//		//.......................................................................................................
       

        //-------------------- IMRT collison where body force has NO higher-order terms -------------//
		//..............carry out relaxation process...............................................
		m1 = m1 + rlx_setA*((-30*Den+19*Den*(ux*ux+uy*uy+uz*uz)/porosity + 57*pressure*porosity) - m1);
		m2 = m2 + rlx_setA*((12*Den - 5.5*Den*(ux*ux+uy*uy+uz*uz)/porosity-27*pressure*porosity) - m2);
        jx = jx + Fx;
		m4 = m4 + rlx_setB*((-0.6666666666666666*ux*Den) - m4)
                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
        jy = jy + Fy;
		m6 = m6 + rlx_setB*((-0.6666666666666666*uy*Den) - m6)
                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
        jz = jz + Fz;
		m8 = m8 + rlx_setB*((-0.6666666666666666*uz*Den) - m8)
                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
		m9 = m9 + rlx_setA*((Den*(2*ux*ux-uy*uy-uz*uz)/porosity) - m9);
		m10 = m10 + rlx_setA*(-0.5*Den*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10);
		m11 = m11 + rlx_setA*((Den*(uy*uy-uz*uz)/porosity) - m11);
		m12 = m12 + rlx_setA*(-0.5*(Den*(uy*uy-uz*uz)/porosity)- m12);
		m13 = m13 + rlx_setA*((Den*ux*uy/porosity) - m13);
		m14 = m14 + rlx_setA*((Den*uy*uz/porosity) - m14);
		m15 = m15 + rlx_setA*((Den*ux*uz/porosity) - m15);
		m16 = m16 + rlx_setB*( - m16);
		m17 = m17 + rlx_setB*( - m17);
		m18 = m18 + rlx_setB*( - m18);
		//.......................................................................................................


		//.................inverse transformation......................................................
		// q=0
		fq = mrt_V1*Den-mrt_V2*m1+mrt_V3*m2;
		dist[n] = fq;

		// q = 1
		fq = mrt_V1*Den-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10);
		nread = neighborList[n+Np];
		dist[nread] = fq;

		// q=2
		fq = mrt_V1*Den-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10);
		nread = neighborList[n];
		dist[nread] = fq;

		// q = 3
		fq = mrt_V1*Den-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
		nread = neighborList[n+3*Np];
		dist[nread] = fq;

		// q = 4
		fq = mrt_V1*Den-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
		nread = neighborList[n+2*Np];
		dist[nread] = fq;

		// q = 5
		fq = mrt_V1*Den-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
		nread = neighborList[n+5*Np];
		dist[nread] = fq;

		// q = 6
		fq = mrt_V1*Den-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
		nread = neighborList[n+4*Np];
		dist[nread] = fq;

		// q = 7
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m16-m17);
		nread = neighborList[n+7*Np];
		dist[nread] = fq;

		// q = 8
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m17-m16);
		nread = neighborList[n+6*Np];
		dist[nread] = fq;

		// q = 9
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13+0.125*(m16+m17);
		nread = neighborList[n+9*Np];
		dist[nread] = fq;

		// q = 10
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13-0.125*(m16+m17);
		nread = neighborList[n+8*Np];
		dist[nread] = fq;

		// q = 11
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m18-m16);
		nread = neighborList[n+11*Np];
		dist[nread] = fq;

		// q = 12
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m16-m18);
		nread = neighborList[n+10*Np];
		dist[nread]= fq;

		// q = 13
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12-0.25*m15-0.125*(m16+m18);
		nread = neighborList[n+13*Np];
		dist[nread] = fq;

		// q= 14
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12-0.25*m15+0.125*(m16+m18);
		nread = neighborList[n+12*Np];
		dist[nread] = fq;

		// q = 15
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18);
		nread = neighborList[n+15*Np];
		dist[nread] = fq;

		// q = 16
		fq =  mrt_V1*Den+mrt_V9*m1+mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17);
		nread = neighborList[n+14*Np];
		dist[nread] = fq;

		// q = 17
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)-mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18);
		nread = neighborList[n+17*Np];
		dist[nread] = fq;

		// q = 18
		fq = mrt_V1*Den+mrt_V9*m1+mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)-mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18);
		nread = neighborList[n+16*Np];
		dist[nread] = fq;
		//........................................................................

        //Update velocity on device
		Velocity[0*Np+n] = ux;
		Velocity[1*Np+n] = uy;
		Velocity[2*Np+n] = uz;
        //Update pressure on device
        Pressure[n] = pressure;
	}
}


extern "C" void ScaLBL_D3Q19_AAodd_Greyscale_MRT(int *neighborList, double *dist, int start, int finish, int Np, double rlx, double rlx_eff, double Gx, double Gy, double Gz,double *Poros,double *Perm, double *Velocity,double rho0,double *Pressure){

	int n, nread;
	int nr1,nr2,nr3,nr4,nr5,nr6;
	int nr7,nr8,nr9,nr10;
	int nr11,nr12,nr13,nr14;
	double vx,vy,vz,v_mag;
    double ux,uy,uz,u_mag;
    double pressure;//defined for this incompressible model
	// conserved momemnts
	double rho,jx,jy,jz;
	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
    double fq;
	//double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18;
    double GeoFun;//geometric function from Guo's PRE 66, 036304 (2002)
    double porosity;
    double perm;//voxel permeability
    double c0, c1; //Guo's model parameters
    double mu_eff = (1.0/rlx_eff-0.5)/3.0;//kinematic viscosity
    double Fx, Fy, Fz;//The total body force including Brinkman force and user-specified (Gx,Gy,Gz)
    double rlx_setA = rlx;
    double rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);

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
        //........................................................................
        //					READ THE DISTRIBUTIONS
        //		(read from opposite array due to previous swap operation)
        //........................................................................
        // q=0
		fq = dist[n];
		rho = fq;
		m1  = -30.0*fq;
		m2  = 12.0*fq;

		// q=1
		//nread = neighborList[n]; // neighbor 2 
		//fq = dist[nread]; // reading the f1 data into register fq		
		nr1 = neighborList[n]; 
		fq = dist[nr1]; // reading the f1 data into register fq
		rho += fq;
		m1 -= 11.0*fq;
		m2 -= 4.0*fq;
		jx = fq;
		m4 = -4.0*fq;
		m9 = 2.0*fq;
		m10 = -4.0*fq;

		// f2 = dist[10*Np+n];
		//nread = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
		//fq = dist[nread];  // reading the f2 data into register fq
		nr2 = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
		fq = dist[nr2];  // reading the f2 data into register fq
		rho += fq;
		m1 -= 11.0*(fq);
		m2 -= 4.0*(fq);
		jx -= fq;
		m4 += 4.0*(fq);
		m9 += 2.0*(fq);
		m10 -= 4.0*(fq);

		// q=3
		//nread = neighborList[n+2*Np]; // neighbor 4
		//fq = dist[nread];
		nr3 = neighborList[n+2*Np]; // neighbor 4
		fq = dist[nr3];
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
		//nread = neighborList[n+3*Np]; // neighbor 3
		//fq = dist[nread];
		nr4 = neighborList[n+3*Np]; // neighbor 3
		fq = dist[nr4];
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
		//nread = neighborList[n+4*Np];
		//fq = dist[nread];
		nr5 = neighborList[n+4*Np];
		fq = dist[nr5];
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
		//nread = neighborList[n+5*Np];
		//fq = dist[nread];
		nr6 = neighborList[n+5*Np];
		fq = dist[nr6];
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
		//nread = neighborList[n+6*Np];
		//fq = dist[nread];
		nr7 = neighborList[n+6*Np];
		fq = dist[nr7];
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
		//nread = neighborList[n+7*Np];
		//fq = dist[nread];
		nr8 = neighborList[n+7*Np];
		fq = dist[nr8];
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
		//nread = neighborList[n+8*Np];
		//fq = dist[nread];
		nr9 = neighborList[n+8*Np];
		fq = dist[nr9];
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
		//nread = neighborList[n+9*Np];
		//fq = dist[nread];
		nr10 = neighborList[n+9*Np];
		fq = dist[nr10];
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
		//nread = neighborList[n+10*Np];
		//fq = dist[nread];
		nr11 = neighborList[n+10*Np];
		fq = dist[nr11];
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
		//nread = neighborList[n+11*Np];
		//fq = dist[nread];
		nr12 = neighborList[n+11*Np];
		fq = dist[nr12];
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
		//nread = neighborList[n+12*Np];
		//fq = dist[nread];
		nr13 = neighborList[n+12*Np];
		fq = dist[nr13];
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
		//nread = neighborList[n+13*Np];
		//fq = dist[nread];
		nr14 = neighborList[n+13*Np];
		fq = dist[nr14];
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
        //---------------------------------------------------------------------//

        porosity = Poros[n];
        perm = Perm[n];

        c0 = 0.5*(1.0+porosity*0.5*mu_eff/perm);
        if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
        GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
        c1 = porosity*0.5*GeoFun/sqrt(perm);
        if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes

        vx = jx/rho0+0.5*porosity*Gx;
        vy = jy/rho0+0.5*porosity*Gy;
        vz = jz/rho0+0.5*porosity*Gz;
        v_mag=sqrt(vx*vx+vy*vy+vz*vz);
        ux = vx/(c0+sqrt(c0*c0+c1*v_mag));
        uy = vy/(c0+sqrt(c0*c0+c1*v_mag));
        uz = vz/(c0+sqrt(c0*c0+c1*v_mag));
        u_mag=sqrt(ux*ux+uy*uy+uz*uz);

        //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
        Fx = rho0*(-porosity*mu_eff/perm*ux - porosity*GeoFun/sqrt(perm)*u_mag*ux + porosity*Gx);
        Fy = rho0*(-porosity*mu_eff/perm*uy - porosity*GeoFun/sqrt(perm)*u_mag*uy + porosity*Gy);
        Fz = rho0*(-porosity*mu_eff/perm*uz - porosity*GeoFun/sqrt(perm)*u_mag*uz + porosity*Gz);
        if (porosity==1.0){
            Fx=rho0*Gx;
            Fy=rho0*Gy;
            Fz=rho0*Gz;
        }

        //Calculate pressure for MRT model
        pressure=rho/3.f;

        //-------------------- MRT collison where body force has NO higher-order terms -------------//
		m1 = m1 + rlx_setA*((19*(ux*ux+uy*uy+uz*uz)*rho0/porosity - 11*rho) - m1);
		m2 = m2 + rlx_setA*((3*rho - 5.5*(ux*ux+uy*uy+uz*uz)*rho0/porosity) - m2);
        jx = jx + Fx;
		m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0)- m4)
                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
        jy = jy + Fy;
		m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0)- m6)
                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
        jz = jz + Fz;
		m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0)- m8)
                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
		m9 = m9 + rlx_setA*(((2*ux*ux-uy*uy-uz*uz)*rho0/porosity) - m9);
		m10 = m10 + rlx_setA*( - m10);
        //m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10);
		m11 = m11 + rlx_setA*(((uy*uy-uz*uz)*rho0/porosity) - m11);
		m12 = m12 + rlx_setA*( - m12);
        //m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12);
		m13 = m13 + rlx_setA*( (ux*uy*rho0/porosity) - m13);
		m14 = m14 + rlx_setA*( (uy*uz*rho0/porosity) - m14);
		m15 = m15 + rlx_setA*( (ux*uz*rho0/porosity) - m15);
		m16 = m16 + rlx_setB*( - m16);
		m17 = m17 + rlx_setB*( - m17);
		m18 = m18 + rlx_setB*( - m18);
        //.......................................................................................................


        //.................inverse transformation......................................................
        // q=0
		fq = mrt_V1*rho-mrt_V2*m1+mrt_V3*m2;
		dist[n] = fq;

		// q = 1
		fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10);
		//nread = neighborList[n+Np];
		dist[nr2] = fq;

		// q=2
		fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10);
		//nread = neighborList[n];
		dist[nr1] = fq;

		// q = 3
		fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
		//nread = neighborList[n+3*Np];
		dist[nr4] = fq;

		// q = 4
		fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
		//nread = neighborList[n+2*Np];
		dist[nr3] = fq;

		// q = 5
		fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
		//nread = neighborList[n+5*Np];
		dist[nr6] = fq;

		// q = 6
		fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
		//nread = neighborList[n+4*Np];
		dist[nr5] = fq;

		// q = 7
		fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)+
				mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m16-m17);
		//nread = neighborList[n+7*Np];
		dist[nr8] = fq;

		// q = 8
		fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
				+mrt_V12*m12+0.25*m13+0.125*(m17-m16);
		//nread = neighborList[n+6*Np];
		dist[nr7] = fq;

		// q = 9
		fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)+
				mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13+0.125*(m16+m17);
		//nread = neighborList[n+9*Np];
		dist[nr10] = fq;

		// q = 10
		fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)+
				mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13-0.125*(m16+m17);
		//nread = neighborList[n+8*Np];
		dist[nr9] = fq;

		// q = 11
		fq = mrt_V1*rho+mrt_V9*m1
				+mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)
				+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
				-mrt_V12*m12+0.25*m15+0.125*(m18-m16);
		//nread = neighborList[n+11*Np];
		dist[nr12] = fq;

		// q = 12
		fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)+
				mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m16-m18);
		//nread = neighborList[n+10*Np];
		dist[nr11]= fq;

		// q = 13
		fq = mrt_V1*rho+mrt_V9*m1
				+mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)
				+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
				-mrt_V12*m12-0.25*m15-0.125*(m16+m18);
		//nread = neighborList[n+13*Np];
		dist[nr14] = fq;

		// q= 14
		fq = mrt_V1*rho+mrt_V9*m1
				+mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)
				+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
				-mrt_V12*m12-0.25*m15+0.125*(m16+m18);
		//nread = neighborList[n+12*Np];
		dist[nr13] = fq;


		// q = 15
		fq = mrt_V1*rho+mrt_V9*m1
				+mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)
				-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18);
		nread = neighborList[n+15*Np];
		dist[nread] = fq;

		// q = 16
		fq =  mrt_V1*rho+mrt_V9*m1
				+mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)
				-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17);
		nread = neighborList[n+14*Np];
		dist[nread] = fq;


		// q = 17
		fq = mrt_V1*rho+mrt_V9*m1
				+mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)
				-mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18);
		nread = neighborList[n+17*Np];
		dist[nread] = fq;

		// q = 18
		fq = mrt_V1*rho+mrt_V9*m1
				+mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)
				-mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18);
		nread = neighborList[n+16*Np];
		dist[nread] = fq;
        //........................................................................

        //Update velocity on device
        Velocity[0*Np+n] = ux;
        Velocity[1*Np+n] = uy;
        Velocity[2*Np+n] = uz;
        //Update pressure on device
        Pressure[n] = pressure;
	}
}

extern "C" void ScaLBL_D3Q19_AAeven_Greyscale_MRT(double *dist, int start, int finish, int Np, double rlx, double rlx_eff, double Gx, double Gy, double Gz,double *Poros,double *Perm, double *Velocity,double rho0,double *Pressure){

	int n;
	double vx,vy,vz,v_mag;
    double ux,uy,uz,u_mag;
    double pressure;//defined for this incompressible model
	// conserved momemnts
	double rho,jx,jy,jz;
	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
    double fq;
	//double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18;
    double GeoFun;//geometric function from Guo's PRE 66, 036304 (2002)
    double porosity;
    double perm;//voxel permeability
    double c0, c1; //Guo's model parameters
    double mu_eff = (1.0/rlx_eff-0.5)/3.0;//kinematic viscosity
    double Fx, Fy, Fz;//The total body force including Brinkman force and user-specified (Gx,Gy,Gz)
    double rlx_setA = rlx;
    double rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);

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
        //........................................................................
        //					READ THE DISTRIBUTIONS
        //		(read from opposite array due to previous swap operation)
        //........................................................................
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
        //---------------------------------------------------------------------//

        porosity = Poros[n];
        perm = Perm[n];

        c0 = 0.5*(1.0+porosity*0.5*mu_eff/perm);
        if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
        GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
        c1 = porosity*0.5*GeoFun/sqrt(perm);
        if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes

        vx = jx/rho0+0.5*porosity*Gx;
        vy = jy/rho0+0.5*porosity*Gy;
        vz = jz/rho0+0.5*porosity*Gz;
        v_mag=sqrt(vx*vx+vy*vy+vz*vz);
        ux = vx/(c0+sqrt(c0*c0+c1*v_mag));
        uy = vy/(c0+sqrt(c0*c0+c1*v_mag));
        uz = vz/(c0+sqrt(c0*c0+c1*v_mag));
        u_mag=sqrt(ux*ux+uy*uy+uz*uz);

        //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
        Fx = rho0*(-porosity*mu_eff/perm*ux - porosity*GeoFun/sqrt(perm)*u_mag*ux + porosity*Gx);
        Fy = rho0*(-porosity*mu_eff/perm*uy - porosity*GeoFun/sqrt(perm)*u_mag*uy + porosity*Gy);
        Fz = rho0*(-porosity*mu_eff/perm*uz - porosity*GeoFun/sqrt(perm)*u_mag*uz + porosity*Gz);
        if (porosity==1.0){
            Fx=rho0*Gx;
            Fy=rho0*Gy;
            Fz=rho0*Gz;
        }

        //Calculate pressure for Incompressible-MRT model
        pressure=rho/3.f;

        //-------------------- IMRT collison where body force has NO higher-order terms -------------//
		m1 = m1 + rlx_setA*((19*(ux*ux+uy*uy+uz*uz)*rho0/porosity - 11*rho) - m1);
		m2 = m2 + rlx_setA*((3*rho - 5.5*(ux*ux+uy*uy+uz*uz)*rho0/porosity) - m2);
        jx = jx + Fx;
		m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0)- m4)
                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
        jy = jy + Fy;
		m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0)- m6)
                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
        jz = jz + Fz;
		m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0)- m8)
                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
		m9 = m9 + rlx_setA*(((2*ux*ux-uy*uy-uz*uz)*rho0/porosity) - m9);
		m10 = m10 + rlx_setA*( - m10);
        //m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10);
		m11 = m11 + rlx_setA*(((uy*uy-uz*uz)*rho0/porosity) - m11);
		m12 = m12 + rlx_setA*( - m12);
        //m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12);
		m13 = m13 + rlx_setA*( (ux*uy*rho0/porosity) - m13);
		m14 = m14 + rlx_setA*( (uy*uz*rho0/porosity) - m14);
		m15 = m15 + rlx_setA*( (ux*uz*rho0/porosity) - m15);
		m16 = m16 + rlx_setB*( - m16);
		m17 = m17 + rlx_setB*( - m17);
		m18 = m18 + rlx_setB*( - m18);
        //.......................................................................................................

        //.................inverse transformation......................................................
        // q=0
		fq = mrt_V1*rho-mrt_V2*m1+mrt_V3*m2;
		dist[n] = fq;

		// q = 1
		fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10);
		dist[1*Np+n] = fq;

		// q=2
		fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10);
		dist[2*Np+n] = fq;

		// q = 3
		fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
		dist[3*Np+n] = fq;

		// q = 4
		fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
		dist[4*Np+n] = fq;

		// q = 5
		fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
		dist[5*Np+n] = fq;

		// q = 6
		fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
		dist[6*Np+n] = fq;

		// q = 7
		fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)+
				mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m16-m17);
		dist[7*Np+n] = fq;


		// q = 8
		fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
				+mrt_V12*m12+0.25*m13+0.125*(m17-m16);
		dist[8*Np+n] = fq;

		// q = 9
		fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)+
				mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13+0.125*(m16+m17);
		dist[9*Np+n] = fq;

		// q = 10
		fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)+
				mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13-0.125*(m16+m17);
		dist[10*Np+n] = fq;


		// q = 11
		fq = mrt_V1*rho+mrt_V9*m1
				+mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)
				+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
				-mrt_V12*m12+0.25*m15+0.125*(m18-m16);
		dist[11*Np+n] = fq;

		// q = 12
		fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)+
				mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m16-m18);
		dist[12*Np+n] = fq;

		// q = 13
		fq = mrt_V1*rho+mrt_V9*m1
				+mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)
				+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
				-mrt_V12*m12-0.25*m15-0.125*(m16+m18);
		dist[13*Np+n] = fq;

		// q= 14
		fq = mrt_V1*rho+mrt_V9*m1
				+mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)
				+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
				-mrt_V12*m12-0.25*m15+0.125*(m16+m18);

		dist[14*Np+n] = fq;

		// q = 15
		fq = mrt_V1*rho+mrt_V9*m1
				+mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)
				-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18);
		dist[15*Np+n] = fq;

		// q = 16
		fq =  mrt_V1*rho+mrt_V9*m1
				+mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)
				-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17);
		dist[16*Np+n] = fq;


		// q = 17
		fq = mrt_V1*rho+mrt_V9*m1
				+mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)
				-mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18);
		dist[17*Np+n] = fq;

		// q = 18
		fq = mrt_V1*rho+mrt_V9*m1
				+mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)
				-mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18);
		dist[18*Np+n] = fq;
        //........................................................................

        //Update velocity on device
        Velocity[0*Np+n] = ux;
        Velocity[1*Np+n] = uy;
        Velocity[2*Np+n] = uz;
        //Update pressure on device
        Pressure[n] = pressure;
	}
}


extern "C" void ScaLBL_D3Q19_GreyIMRT_Init(double *dist, int Np, double Den)
{
	int n;
	for (n=0; n<Np; n++){
		dist[n] = Den - 0.6666666666666667;
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

