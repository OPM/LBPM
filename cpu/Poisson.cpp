
extern "C" void ScaLBL_D3Q7_AAodd_Poisson(int *neighborList, double *dist, double *Den_charge, double *Psi, double *ElectricField, double tau, double epsilon_LB,double gamma,
        int start, int finish, int Np){
	int n;
	double psi;//electric potential
    double Ex,Ey,Ez;//electrical field
    double rho_e;//local charge density
	double f0,f1,f2,f3,f4,f5,f6;
	int nr1,nr2,nr3,nr4,nr5,nr6;
    double rlx=1.0/tau;

	for (n=start; n<finish; n++){

        //Load data
        rho_e = Den_charge[n];
        rho_e = gamma*rho_e/epsilon_LB;
		
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
		
		psi = f0+f2+f1+f4+f3+f6+f5;
		Ex = (f1-f2)*rlx*4.5;//NOTE the unit of electric field here is V/lu
		Ey = (f3-f4)*rlx*4.5;
		Ez = (f5-f6)*rlx*4.5;
        ElectricField[n+0*Np] = Ex;
        ElectricField[n+1*Np] = Ey;
        ElectricField[n+2*Np] = Ez;
        Psi[n] = psi;

		// q = 0
		dist[n] = f0*(1.0-rlx) + 0.3333333333333333*(rlx*psi+rho_e);

		// q = 1
		dist[nr2] = f1*(1.0-rlx) + 0.1111111111111111*(rlx*psi+rho_e);

		// q = 2
		dist[nr1] = f2*(1.0-rlx) + 0.1111111111111111*(rlx*psi+rho_e);

		// q = 3
		dist[nr4] = f3*(1.0-rlx) + 0.1111111111111111*(rlx*psi+rho_e);

		// q = 4
		dist[nr3] = f4*(1.0-rlx) + 0.1111111111111111*(rlx*psi+rho_e);

		// q = 5
		dist[nr6] = f5*(1.0-rlx) + 0.1111111111111111*(rlx*psi+rho_e);

		// q = 6
		dist[nr5] = f6*(1.0-rlx) + 0.1111111111111111*(rlx*psi+rho_e);
		//........................................................................
	}
}

extern "C" void ScaLBL_D3Q7_AAeven_Poisson(double *dist, double *Den_charge, double *Psi, double *ElectricField, double tau, double epsilon_LB,double gamma,
        int start, int finish, int Np){
	int n;
	double psi;//electric potential
    double Ex,Ey,Ez;//electrical field
    double rho_e;//local charge density
	double f0,f1,f2,f3,f4,f5,f6;
    double rlx=1.0/tau;

	for (n=start; n<finish; n++){

        //Load data
        rho_e = Den_charge[n];
        rho_e = gamma*rho_e/epsilon_LB;

		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f5 = dist[6*Np+n];
		f6 = dist[5*Np+n];


		psi = f0+f2+f1+f4+f3+f6+f5;
		Ex = (f1-f2)*rlx*4.5;//NOTE the unit of electric field here is V/lu
		Ey = (f3-f4)*rlx*4.5;
		Ez = (f5-f6)*rlx*4.5;
        ElectricField[n+0*Np] = Ex;
        ElectricField[n+1*Np] = Ey;
        ElectricField[n+2*Np] = Ez;
        Psi[n] = psi;

		// q = 0
		dist[n] = f0*(1.0-rlx) + 0.3333333333333333*(rlx*psi+rho_e);

		// q = 1
		dist[1*Np+n] = f1*(1.0-rlx) + 0.1111111111111111*(rlx*psi+rho_e); 

		// q = 2
		dist[2*Np+n] = f2*(1.0-rlx) + 0.1111111111111111*(rlx*psi+rho_e); 

		// q = 3
		dist[3*Np+n] = f3*(1.0-rlx) + 0.1111111111111111*(rlx*psi+rho_e);

		// q = 4
		dist[4*Np+n] = f4*(1.0-rlx) + 0.1111111111111111*(rlx*psi+rho_e); 

		// q = 5
		dist[5*Np+n] = f5*(1.0-rlx) + 0.1111111111111111*(rlx*psi+rho_e); 

		// q = 6
		dist[6*Np+n] = f6*(1.0-rlx) + 0.1111111111111111*(rlx*psi+rho_e); 
		//........................................................................
	}
}

extern "C" void ScaLBL_D3Q7_Poisson_Init(double *dist, int Np)
{
	int n;
	for (n=0; n<Np; n++){
		dist[0*Np+n] = 0.3333333333333333;
		dist[1*Np+n] = 0.1111111111111111;		
		dist[2*Np+n] = 0.1111111111111111;	
		dist[3*Np+n] = 0.1111111111111111;	
		dist[4*Np+n] = 0.1111111111111111;	
		dist[5*Np+n] = 0.1111111111111111;	
		dist[6*Np+n] = 0.1111111111111111;	
	}
}
