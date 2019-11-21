#include <math.h>

extern "C" void ScaLBL_D3Q19_AAeven_Greyscale(double *dist, int start, int finish, int Np, double rlx, double Fx, double Fy, double Fz,
                                              double *Poros,double *Perm, double *Velocity){
	int n;
	// conserved momemnts
	double rho,vx,vy,vz,v_mag;
    double ux,uy,uz,u_mag;
    //double uu;
	// non-conserved moments
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18;
    double GeoFun;//geometric function from Guo's PRE 66, 036304 (2002)
    double porosity;
    double perm;//voxel permeability
    double c0, c1; //Guo's model parameters
    double mu = (1.0/rlx-0.5)/3.0;//kinematic viscosity

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

        c0 = 0.5*(1.0+porosity*0.5*mu/perm);
        if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
        GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
        c1 = porosity*0.5*GeoFun/sqrt(perm);
        if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes

		rho = f0+f2+f1+f4+f3+f6+f5+f8+f7+f10+f9+f12+f11+f14+f13+f16+f15+f18+f17;
		vx = (f1-f2+f7-f8+f9-f10+f11-f12+f13-f14)/rho+0.5*porosity*Fx;
		vy = (f3-f4+f7-f8-f9+f10+f15-f16+f17-f18)/rho+0.5*porosity*Fy;
		vz = (f5-f6+f11-f12-f13+f14+f15-f16-f17+f18)/rho+0.5*porosity*Fz;
        v_mag=sqrt(vx*vx+vy*vy+vz*vz);
        ux = vx/(c0+sqrt(c0*c0+c1*v_mag));
        uy = vy/(c0+sqrt(c0*c0+c1*v_mag));
        uz = vz/(c0+sqrt(c0*c0+c1*v_mag));
        u_mag=sqrt(ux*ux+uy*uy+uz*uz);
		//uu = 1.5*(ux*ux+uy*uy+uz*uz);

        //Update the body force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
        double Fx_tmp=Fx; //Fx_tmp stores user-specified body force
        double Fy_tmp=Fy;
        double Fz_tmp=Fz;
        Fx = -porosity*mu/perm*ux - porosity*GeoFun/sqrt(perm)*u_mag*ux + porosity*Fx;
        Fy = -porosity*mu/perm*uy - porosity*GeoFun/sqrt(perm)*u_mag*uy + porosity*Fy;
        Fz = -porosity*mu/perm*uz - porosity*GeoFun/sqrt(perm)*u_mag*uz + porosity*Fz;
        if (porosity==1.0){
            Fx=Fx_tmp;
            Fy=Fy_tmp;
            Fz=Fz_tmp;
        }

		// q=0
		dist[n] = f0*(1.0-rlx)+ rlx*0.3333333333333333*rho*(1. - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                  + 0.3333333333333333*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(0. - (3.*uy)/porosity) + Fz*(0. - (3.*uz)/porosity));

		// q = 1
		dist[1*Np+n] = f1*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 + 3.*ux + (4.5*ux*ux)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
            +0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(3. + (6.*ux)/porosity) + Fy*(0. - (3.*uy)/porosity) + Fz*(0. - (3.*uz)/porosity));

		// q=2
		dist[2*Np+n] = f2*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 - 3.*ux + (4.5*ux*ux)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
            +0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(-3. + (6.*ux)/porosity) + Fy*(0. - (3.*uy)/porosity) + Fz*(0. - (3.*uz)/porosity));

		// q = 3
		dist[3*Np+n] = f3*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 + 3.*uy + (4.5*uy*uy)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(3. + (6.*uy)/porosity) + Fz*(0. - (3.*uz)/porosity));

		// q = 4
		dist[4*Np+n] = f4*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 - 3.*uy + (4.5*uy*uy)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(-3. + (6.*uy)/porosity) + Fz*(0. - (3.*uz)/porosity));

		// q = 5
		dist[5*Np+n] = f5*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 + 3.*uz + (4.5*uz*uz)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(0. - (3.*uy)/porosity) + Fz*(3. + (6.*uz)/porosity));

		// q = 6
		dist[6*Np+n] = f6*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 - 3.*uz + (4.5*uz*uz)/porosity - (1.5*(ux*ux+ uy*uy + uz*uz))/porosity)
				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(0. - (3.*uy)/porosity) + Fz*(-3. + (6.*uz)/porosity));

		// q = 7
		dist[7*Np+n] = f7*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux + uy) + (4.5*(ux + uy)*(ux + uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(3. - (3.*ux)/porosity + (9.*(ux + uy))/porosity) + Fy*(3. - (3.*uy)/porosity + (9.*(ux + uy))/porosity) + 
  Fz*(0. - (3.*uz)/porosity));

		// q = 8
		dist[8*Np+n] = f8*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux - uy) + (4.5*(-ux - uy)*(-ux - uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(-3. - (3.*ux)/porosity - (9.*(-ux - uy))/porosity) + Fy*(-3. - (9.*(-ux - uy))/porosity - (3.*uy)/porosity) + 
  Fz*(0. - (3.*uz)/porosity));

		// q = 9
		dist[9*Np+n] = f9*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux - uy) + (4.5*(ux - uy)*(ux - uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(3. - (3.*ux)/porosity + (9.*(ux - uy))/porosity) + Fy*(-3. - (9.*(ux - uy))/porosity - (3.*uy)/porosity) + 
  Fz*(0. - (3.*uz)/porosity));

		// q = 10
		dist[10*Np+n] = f10*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux + uy) + (4.5*(-ux + uy)*(-ux + uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(-3. - (3.*ux)/porosity - (9.*(-ux + uy))/porosity) + Fy*(3. - (3.*uy)/porosity + (9.*(-ux + uy))/porosity) + 
  Fz*(0. - (3.*uz)/porosity));

		// q = 11
		dist[11*Np+n] = f11*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux + uz) + (4.5*(ux + uz)*(ux + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(0. - (3.*uy)/porosity) + Fx*(3. - (3.*ux)/porosity + (9.*(ux + uz))/porosity) + 
  Fz*(3. - (3.*uz)/porosity + (9.*(ux + uz))/porosity));

		// q = 12
		dist[12*Np+n] = f12*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux - uz) + (4.5*(-ux - uz)*(-ux - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(0. - (3.*uy)/porosity) + Fx*(-3. - (3.*ux)/porosity - (9.*(-ux - uz))/porosity) + 
  Fz*(-3. - (9.*(-ux - uz))/porosity - (3.*uz)/porosity));

		// q = 13
		dist[13*Np+n] = f13*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux - uz) + (4.5*(ux - uz)*(ux - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(0. - (3.*uy)/porosity) + Fx*(3. - (3.*ux)/porosity + (9.*(ux - uz))/porosity) + 
  Fz*(-3. - (9.*(ux - uz))/porosity - (3.*uz)/porosity));

		// q= 14
		dist[14*Np+n] = f14*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux + uz) + (4.5*(-ux + uz)*(-ux + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(0. - (3.*uy)/porosity) + Fx*(-3. - (3.*ux)/porosity - (9.*(-ux + uz))/porosity) + 
  Fz*(3. - (3.*uz)/porosity + (9.*(-ux + uz))/porosity));

		// q = 15
		dist[15*Np+n] = f15*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(uy + uz) + (4.5*(uy + uz)*(uy + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(3. - (3.*uy)/porosity + (9.*(uy + uz))/porosity) + 
  Fz*(3. - (3.*uz)/porosity + (9.*(uy + uz))/porosity));

		// q = 16
		dist[16*Np+n] = f16*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-uy - uz) + (4.5*(-uy - uz)*(-uy - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(-3. - (3.*uy)/porosity - (9.*(-uy - uz))/porosity) + 
  Fz*(-3. - (9.*(-uy - uz))/porosity - (3.*uz)/porosity));

		// q = 17
		dist[17*Np+n] = f17*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(uy - uz) + (4.5*(uy - uz)*(uy - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(3. - (3.*uy)/porosity + (9.*(uy - uz))/porosity) + 
  Fz*(-3. - (9.*(uy - uz))/porosity - (3.*uz)/porosity));

		// q = 18
		dist[18*Np+n] = f18*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-uy + uz) + (4.5*(-uy + uz)*(-uy + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(-3. - (3.*uy)/porosity - (9.*(-uy + uz))/porosity) + 
  Fz*(3. - (3.*uz)/porosity + (9.*(-uy + uz))/porosity));

        //Update velocity on device
		Velocity[0*Np+n] = ux;
		Velocity[1*Np+n] = uy;
		Velocity[2*Np+n] = uz;
	}
}

extern "C" void ScaLBL_D3Q19_AAodd_Greyscale(int *neighborList, double *dist, int start, int finish, int Np, double rlx, double Fx, double Fy, double Fz, 
                                             double *Poros,double *Perm, double *Velocity){
	int n;
	// conserved momemnts
	double rho,vx,vy,vz,v_mag;
    double ux,uy,uz,u_mag;
    //double uu;
	// non-conserved moments
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18;
	int nr1,nr2,nr3,nr4,nr5,nr6,nr7,nr8,nr9,nr10,nr11,nr12,nr13,nr14,nr15,nr16,nr17,nr18;
    double GeoFun;//geometric function from Guo's PRE 66, 036304 (2002)
    double porosity;
    double perm;//voxel permeability
    double c0, c1; //Guo's model parameters
    double mu = (1.0/rlx-0.5)/3.0;//kinematic viscosity

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

        c0 = 0.5*(1.0+porosity*0.5*mu/perm);
        if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
        GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
        c1 = porosity*0.5*GeoFun/sqrt(perm);
        if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes

		rho = f0+f2+f1+f4+f3+f6+f5+f8+f7+f10+f9+f12+f11+f14+f13+f16+f15+f18+f17;
		vx = (f1-f2+f7-f8+f9-f10+f11-f12+f13-f14)/rho+0.5*porosity*Fx;
		vy = (f3-f4+f7-f8-f9+f10+f15-f16+f17-f18)/rho+0.5*porosity*Fy;
		vz = (f5-f6+f11-f12-f13+f14+f15-f16-f17+f18)/rho+0.5*porosity*Fz;
        v_mag=sqrt(vx*vx+vy*vy+vz*vz);
        ux = vx/(c0+sqrt(c0*c0+c1*v_mag));
        uy = vy/(c0+sqrt(c0*c0+c1*v_mag));
        uz = vz/(c0+sqrt(c0*c0+c1*v_mag));
        u_mag=sqrt(ux*ux+uy*uy+uz*uz);
		//uu = 1.5*(ux*ux+uy*uy+uz*uz);

        //Update the body force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
        double Fx_tmp=Fx; //Fx_tmp stores user-specified body force
        double Fy_tmp=Fy;
        double Fz_tmp=Fz;
        Fx = -porosity*mu/perm*ux - porosity*GeoFun/sqrt(perm)*u_mag*ux + porosity*Fx;
        Fy = -porosity*mu/perm*uy - porosity*GeoFun/sqrt(perm)*u_mag*uy + porosity*Fy;
        Fz = -porosity*mu/perm*uz - porosity*GeoFun/sqrt(perm)*u_mag*uz + porosity*Fz;
        if (porosity==1.0){
            Fx=Fx_tmp;
            Fy=Fy_tmp;
            Fz=Fz_tmp;
        }

		// q=0
		dist[n] = f0*(1.0-rlx) + rlx*0.3333333333333333*rho*(1. - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                  + 0.3333333333333333*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(0. - (3.*uy)/porosity) + Fz*(0. - (3.*uz)/porosity));

		// q = 1
		dist[nr2] = f1*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 + 3.*ux + (4.5*ux*ux)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
            +0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(3. + (6.*ux)/porosity) + Fy*(0. - (3.*uy)/porosity) + Fz*(0. - (3.*uz)/porosity));

		// q=2
		dist[nr1] = f2*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 - 3.*ux + (4.5*ux*ux)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
            +0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(-3. + (6.*ux)/porosity) + Fy*(0. - (3.*uy)/porosity) + Fz*(0. - (3.*uz)/porosity));

		// q = 3
		dist[nr4] = f3*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 + 3.*uy + (4.5*uy*uy)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(3. + (6.*uy)/porosity) + Fz*(0. - (3.*uz)/porosity));

		// q = 4
		dist[nr3] = f4*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 - 3.*uy + (4.5*uy*uy)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)  
				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(-3. + (6.*uy)/porosity) + Fz*(0. - (3.*uz)/porosity));

		// q = 5
		dist[nr6] = f5*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 + 3.*uz + (4.5*uz*uz)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(0. - (3.*uy)/porosity) + Fz*(3. + (6.*uz)/porosity));

		// q = 6
		dist[nr5] = f6*(1.0-rlx) + rlx*0.05555555555555555*rho*(1 - 3.*uz + (4.5*uz*uz)/porosity - (1.5*(ux*ux+ uy*uy + uz*uz))/porosity) 
				+0.05555555555555555*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(0. - (3.*uy)/porosity) + Fz*(-3. + (6.*uz)/porosity));

		// q = 7
		dist[nr8] = f7*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux + uy) + (4.5*(ux + uy)*(ux + uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(3. - (3.*ux)/porosity + (9.*(ux + uy))/porosity) + Fy*(3. - (3.*uy)/porosity + (9.*(ux + uy))/porosity) + 
  Fz*(0. - (3.*uz)/porosity));

		// q = 8
		dist[nr7] = f8*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux - uy) + (4.5*(-ux - uy)*(-ux - uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(-3. - (3.*ux)/porosity - (9.*(-ux - uy))/porosity) + Fy*(-3. - (9.*(-ux - uy))/porosity - (3.*uy)/porosity) + 
  Fz*(0. - (3.*uz)/porosity));

		// q = 9
		dist[nr10] = f9*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux - uy) + (4.5*(ux - uy)*(ux - uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(3. - (3.*ux)/porosity + (9.*(ux - uy))/porosity) + Fy*(-3. - (9.*(ux - uy))/porosity - (3.*uy)/porosity) + 
  Fz*(0. - (3.*uz)/porosity));

		// q = 10
		dist[nr9] = f10*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux + uy) + (4.5*(-ux + uy)*(-ux + uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(-3. - (3.*ux)/porosity - (9.*(-ux + uy))/porosity) + Fy*(3. - (3.*uy)/porosity + (9.*(-ux + uy))/porosity) + 
  Fz*(0. - (3.*uz)/porosity));

		// q = 11
		dist[nr12] = f11*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux + uz) + (4.5*(ux + uz)*(ux + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(0. - (3.*uy)/porosity) + Fx*(3. - (3.*ux)/porosity + (9.*(ux + uz))/porosity) + 
  Fz*(3. - (3.*uz)/porosity + (9.*(ux + uz))/porosity));

		// q = 12
		dist[nr11] = f12*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux - uz) + (4.5*(-ux - uz)*(-ux - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(0. - (3.*uy)/porosity) + Fx*(-3. - (3.*ux)/porosity - (9.*(-ux - uz))/porosity) + 
  Fz*(-3. - (9.*(-ux - uz))/porosity - (3.*uz)/porosity));

		// q = 13
		dist[nr14] = f13*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(ux - uz) + (4.5*(ux - uz)*(ux - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(0. - (3.*uy)/porosity) + Fx*(3. - (3.*ux)/porosity + (9.*(ux - uz))/porosity) + 
  Fz*(-3. - (9.*(ux - uz))/porosity - (3.*uz)/porosity));

		// q= 14
		dist[nr13] = f14*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-ux + uz) + (4.5*(-ux + uz)*(-ux + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fy*(0. - (3.*uy)/porosity) + Fx*(-3. - (3.*ux)/porosity - (9.*(-ux + uz))/porosity) + 
  Fz*(3. - (3.*uz)/porosity + (9.*(-ux + uz))/porosity));

		// q = 15
		dist[nr16] = f15*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(uy + uz) + (4.5*(uy + uz)*(uy + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(3. - (3.*uy)/porosity + (9.*(uy + uz))/porosity) + 
  Fz*(3. - (3.*uz)/porosity + (9.*(uy + uz))/porosity));

		// q = 16
		dist[nr15] = f16*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-uy - uz) + (4.5*(-uy - uz)*(-uy - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(-3. - (3.*uy)/porosity - (9.*(-uy - uz))/porosity) + 
  Fz*(-3. - (9.*(-uy - uz))/porosity - (3.*uz)/porosity));

		// q = 17
		dist[nr18] = f17*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(uy - uz) + (4.5*(uy - uz)*(uy - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(3. - (3.*uy)/porosity + (9.*(uy - uz))/porosity) + 
  Fz*(-3. - (9.*(uy - uz))/porosity - (3.*uz)/porosity));

		// q = 18
		dist[nr17] = f18*(1.0-rlx) + rlx*0.027777777777777776*rho*(1 + 3.*(-uy + uz) + (4.5*(-uy + uz)*(-uy + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
				+0.027777777777777776*rho*(1. - 0.5*rlx)*(Fx*(0. - (3.*ux)/porosity) + Fy*(-3. - (3.*uy)/porosity - (9.*(-uy + uz))/porosity) + 
  Fz*(3. - (3.*uz)/porosity + (9.*(-uy + uz))/porosity));

        //Update velocity on device
		Velocity[0*Np+n] = ux;
		Velocity[1*Np+n] = uy;
		Velocity[2*Np+n] = uz;
	}
}
