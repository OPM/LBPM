extern "C" void ScaLBL_D3Q7_AAeven_Ion(double *dist, double *Velocity, int start, int finish, int Np, double rlx, double Fx, double Fy, double Fz){
	int n;
	// conserved momemnts
	double rho,ux,uy,uz,uu;
	// non-conserved moments
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18;

	for (int n=start; n<finish; n++){
		// q=0
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f5 = dist[6*Np+n];
		f6 = dist[5*Np+n];

		rho = f0+f2+f1+f4+f3+f6;
		ux = Velocity[n];
		uy = Velocity[n+Np];
		uz = Velocity[n+2*Np];
		uu = 1.5*(ux*ux+uy*uy+uz*uz);

		// q=0
		dist[n] = f0*(1.0-rlx)+rlx*0.3333333333333333*(1.0-uu);

		// q = 1
		dist[1*Np+n] = f1*(1.0-rlx) + rlx*0.05555555555555555*(rho + 3.0*ux + 4.5*ux*ux - uu) + 0.16666666*Fx;

		// q=2
		dist[2*Np+n] = f2*(1.0-rlx) + rlx*0.05555555555555555*(rho - 3.0*ux + 4.5*ux*ux - uu)-  0.16666666*Fx;

		// q = 3
		dist[3*Np+n] = f3*(1.0-rlx) +
				rlx*0.05555555555555555*(rho + 3.0*uy + 4.5*uy*uy - uu) + 0.16666666*Fy;

		// q = 4
		dist[4*Np+n] = f4*(1.0-rlx) + 
				rlx*0.05555555555555555*(rho - 3.0*uy + 4.5*uy*uy - uu)- 0.16666666*Fy;

		// q = 5
		dist[5*Np+n] = f5*(1.0-rlx) + 
				rlx*0.05555555555555555*(rho + 3.0*uz + 4.5*uz*uz - uu) + 0.16666666*Fz;

		// q = 6
		dist[6*Np+n] = f6*(1.0-rlx) + 
				rlx*0.05555555555555555*(rho - 3.0*uz + 4.5*uz*uz - uu) - 0.16666666*Fz;


		//........................................................................
	}
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion(int *neighborList, double *dist, double *Velocity, int start, int finish, int Np, double rlx, double Fx, double Fy, double Fz){
	int n;
	// conserved momemnts
	double rho,ux,uy,uz,uu;
	// non-conserved moments
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18;
	int nr1,nr2,nr3,nr4,nr5,nr6,nr7,nr8,nr9,nr10,nr11,nr12,nr13,nr14,nr15,nr16,nr17,nr18;

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
		
		rho = f0+f2+f1+f4+f3+f6;
		ux = Velocity[n];
		uy = Velocity[n+Np];
		uz = Velocity[n+2*Np];
		uu = 1.5*(ux*ux+uy*uy+uz*uz);

		// q=0
		dist[n] = f0*(1.0-rlx)+rlx*0.3333333333333333*(1.0-uu);

		// q = 1
		dist[nr2] = f1*(1.0-rlx) + rlx*0.05555555555555555*(rho + 3.0*ux + 4.5*ux*ux - uu) + 0.16666666*Fx;

		// q=2
		dist[nr1] = f2*(1.0-rlx) + rlx*0.05555555555555555*(rho - 3.0*ux + 4.5*ux*ux - uu)-  0.16666666*Fx;

		// q = 3
		dist[nr4] = f3*(1.0-rlx) +
				rlx*0.05555555555555555*(rho + 3.0*uy + 4.5*uy*uy - uu) + 0.16666666*Fy;

		// q = 4
		dist[nr3] = f4*(1.0-rlx) + 
				rlx*0.05555555555555555*(rho - 3.0*uy + 4.5*uy*uy - uu)- 0.16666666*Fy;

		// q = 5
		dist[nr6] = f5*(1.0-rlx) + 
				rlx*0.05555555555555555*(rho + 3.0*uz + 4.5*uz*uz - uu) + 0.16666666*Fz;

		// q = 6
		dist[nr5] = f6*(1.0-rlx) + 
				rlx*0.05555555555555555*(rho - 3.0*uz + 4.5*uz*uz - uu) - 0.16666666*Fz;

	
	}
}