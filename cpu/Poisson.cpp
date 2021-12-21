
extern "C" void
ScaLBL_D3Q7_AAodd_Poisson_ElectricPotential(int *neighborList, int *Map,
                                            double *dist, double *Psi,
                                            int start, int finish, int Np) {
    int n;
    double psi; //electric potential
    double fq;
    int nread;
    int idx;

    for (n = start; n < finish; n++) {

        // q=0
        fq = dist[n];
        psi = fq;

        // q=1
        nread = neighborList[n];
        fq = dist[nread];
        psi += fq;

        // q=2
        nread = neighborList[n + Np];
        fq = dist[nread];
        psi += fq;

        // q=3
        nread = neighborList[n + 2 * Np];
        fq = dist[nread];
        psi += fq;

        // q = 4
        nread = neighborList[n + 3 * Np];
        fq = dist[nread];
        psi += fq;

        // q=5
        nread = neighborList[n + 4 * Np];
        fq = dist[nread];
        psi += fq;

        // q = 6
        nread = neighborList[n + 5 * Np];
        fq = dist[nread];
        psi += fq;

        idx = Map[n];
        Psi[idx] = psi;
    }
}

extern "C" void ScaLBL_D3Q7_AAeven_Poisson_ElectricPotential(
    int *Map, double *dist, double *Psi, int start, int finish, int Np) {
    int n;
    double psi; //electric potential
    double fq;
    int idx;

    for (n = start; n < finish; n++) {

        // q=0
        fq = dist[n];
        psi = fq;

        // q=1
        fq = dist[2 * Np + n];
        psi += fq;

        // q=2
        fq = dist[1 * Np + n];
        psi += fq;

        // q=3
        fq = dist[4 * Np + n];
        psi += fq;

        // q=4
        fq = dist[3 * Np + n];
        psi += fq;

        // q=5
        fq = dist[6 * Np + n];
        psi += fq;

        // q=6
        fq = dist[5 * Np + n];
        psi += fq;

        idx = Map[n];
        Psi[idx] = psi;
    }
}

extern "C" void ScaLBL_D3Q7_AAodd_Poisson(int *neighborList, int *Map,
                                          double *dist, double *Den_charge,
                                          double *Psi, double *ElectricField,
                                          double tau, double epsilon_LB,
                                          int start, int finish, int Np) {
    int n;
    double psi;        //electric potential
    double Ex, Ey, Ez; //electric field
    double rho_e;      //local charge density
    double f0, f1, f2, f3, f4, f5, f6;
    int nr1, nr2, nr3, nr4, nr5, nr6;
    double rlx = 1.0 / tau;
    int idx;

    for (n = start; n < finish; n++) {

        //Load data
        rho_e = Den_charge[n];
        rho_e = rho_e / epsilon_LB;
        idx = Map[n];
        psi = Psi[idx];

        // q=0
        f0 = dist[n];
        // q=1
        nr1 = neighborList[n]; // neighbor 2 ( > 10Np => odd part of dist)
        f1 = dist[nr1];        // reading the f1 data into register fq

        nr2 = neighborList[n + Np]; // neighbor 1 ( < 10Np => even part of dist)
        f2 = dist[nr2];             // reading the f2 data into register fq

        // q=3
        nr3 = neighborList[n + 2 * Np]; // neighbor 4
        f3 = dist[nr3];

        // q = 4
        nr4 = neighborList[n + 3 * Np]; // neighbor 3
        f4 = dist[nr4];

        // q=5
        nr5 = neighborList[n + 4 * Np];
        f5 = dist[nr5];

        // q = 6
        nr6 = neighborList[n + 5 * Np];
        f6 = dist[nr6];

        Ex = (f1 - f2) * rlx *
             4.0; //NOTE the unit of electric field here is V/lu
        Ey = (f3 - f4) * rlx *
             4.0; //factor 4.0 is D3Q7 lattice squared speed of sound
        Ez = (f5 - f6) * rlx * 4.0;
        ElectricField[n + 0 * Np] = Ex;
        ElectricField[n + 1 * Np] = Ey;
        ElectricField[n + 2 * Np] = Ez;

        // q = 0
        dist[n] = f0 * (1.0 - rlx) + 0.25 * (rlx * psi + rho_e);

        // q = 1
        dist[nr2] = f1 * (1.0 - rlx) + 0.125 * (rlx * psi + rho_e);

        // q = 2
        dist[nr1] = f2 * (1.0 - rlx) + 0.125 * (rlx * psi + rho_e);

        // q = 3
        dist[nr4] = f3 * (1.0 - rlx) + 0.125 * (rlx * psi + rho_e);

        // q = 4
        dist[nr3] = f4 * (1.0 - rlx) + 0.125 * (rlx * psi + rho_e);

        // q = 5
        dist[nr6] = f5 * (1.0 - rlx) + 0.125 * (rlx * psi + rho_e);

        // q = 6
        dist[nr5] = f6 * (1.0 - rlx) + 0.125 * (rlx * psi + rho_e);
        //........................................................................
    }
}

extern "C" void ScaLBL_D3Q7_AAeven_Poisson(int *Map, double *dist,
                                           double *Den_charge, double *Psi,
                                           double *ElectricField, double tau,
                                           double epsilon_LB, int start,
                                           int finish, int Np) {
    int n;
    double psi;        //electric potential
    double Ex, Ey, Ez; //electric field
    double rho_e;      //local charge density
    double f0, f1, f2, f3, f4, f5, f6;
    double rlx = 1.0 / tau;
    int idx;

    for (n = start; n < finish; n++) {

        //Load data
        rho_e = Den_charge[n];
        rho_e = rho_e / epsilon_LB;
        idx = Map[n];
        psi = Psi[idx];

        f0 = dist[n];
        f1 = dist[2 * Np + n];
        f2 = dist[1 * Np + n];
        f3 = dist[4 * Np + n];
        f4 = dist[3 * Np + n];
        f5 = dist[6 * Np + n];
        f6 = dist[5 * Np + n];

        Ex = (f1 - f2) * rlx *
             4.0; //NOTE the unit of electric field here is V/lu
        Ey = (f3 - f4) * rlx *
             4.0; //factor 4.0 is D3Q7 lattice squared speed of sound
        Ez = (f5 - f6) * rlx * 4.0;
        ElectricField[n + 0 * Np] = Ex;
        ElectricField[n + 1 * Np] = Ey;
        ElectricField[n + 2 * Np] = Ez;

        // q = 0
        dist[n] = f0 * (1.0 - rlx) + 0.25 * (rlx * psi + rho_e);

        // q = 1
        dist[1 * Np + n] = f1 * (1.0 - rlx) + 0.125 * (rlx * psi + rho_e);

        // q = 2
        dist[2 * Np + n] = f2 * (1.0 - rlx) + 0.125 * (rlx * psi + rho_e);

        // q = 3
        dist[3 * Np + n] = f3 * (1.0 - rlx) + 0.125 * (rlx * psi + rho_e);

        // q = 4
        dist[4 * Np + n] = f4 * (1.0 - rlx) + 0.125 * (rlx * psi + rho_e);

        // q = 5
        dist[5 * Np + n] = f5 * (1.0 - rlx) + 0.125 * (rlx * psi + rho_e);

        // q = 6
        dist[6 * Np + n] = f6 * (1.0 - rlx) + 0.125 * (rlx * psi + rho_e);
        //........................................................................
    }
}

extern "C" void ScaLBL_D3Q7_Poisson_Init(int *Map, double *dist, double *Psi,
                                         int start, int finish, int Np) {
    int n;
    int ijk;
    for (n = start; n < finish; n++) {
        ijk = Map[n];
        dist[0 * Np + n] = 0.25 * Psi[ijk];
        dist[1 * Np + n] = 0.125 * Psi[ijk];
        dist[2 * Np + n] = 0.125 * Psi[ijk];
        dist[3 * Np + n] = 0.125 * Psi[ijk];
        dist[4 * Np + n] = 0.125 * Psi[ijk];
        dist[5 * Np + n] = 0.125 * Psi[ijk];
        dist[6 * Np + n] = 0.125 * Psi[ijk];
    }
}

extern "C" void ScaLBL_D3Q7_PoissonResidualError(
    int *neighborList, int *Map, double *ResidualError, double *Psi,
    double *Den_charge, double epsilon_LB, int strideY, int strideZ, int start,
    int finish) {

    int n, nn, ijk;
    double psi;   //electric potential
    double rho_e; //local charge density
    // neighbors of electric potential psi
    double m1, m2, m4, m6, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18;
    double m3, m5, m7;
    double psi_Laplacian;
    double residual_error;

    for (n = start; n < finish; n++) {

        //Load data
        rho_e = Den_charge[n];
        ijk = Map[n];
        psi = Psi[ijk];

        //					COMPUTE THE COLOR GRADIENT
        //........................................................................
        //.................Read Phase Indicator Values............................
        //........................................................................
        nn = ijk - 1; // neighbor index (get convention)
        m1 = Psi[nn]; // get neighbor for phi - 1
        //........................................................................
        nn = ijk + 1; // neighbor index (get convention)
        m2 = Psi[nn]; // get neighbor for phi - 2
        //........................................................................
        nn = ijk - strideY; // neighbor index (get convention)
        m3 = Psi[nn];       // get neighbor for phi - 3
        //........................................................................
        nn = ijk + strideY; // neighbor index (get convention)
        m4 = Psi[nn];       // get neighbor for phi - 4
        //........................................................................
        nn = ijk - strideZ; // neighbor index (get convention)
        m5 = Psi[nn];       // get neighbor for phi - 5
        //........................................................................
        nn = ijk + strideZ; // neighbor index (get convention)
        m6 = Psi[nn];       // get neighbor for phi - 6
        //........................................................................
        nn = ijk - strideY - 1; // neighbor index (get convention)
        m7 = Psi[nn];           // get neighbor for phi - 7
        //........................................................................
        nn = ijk + strideY + 1; // neighbor index (get convention)
        m8 = Psi[nn];           // get neighbor for phi - 8
        //........................................................................
        nn = ijk + strideY - 1; // neighbor index (get convention)
        m9 = Psi[nn];           // get neighbor for phi - 9
        //........................................................................
        nn = ijk - strideY + 1; // neighbor index (get convention)
        m10 = Psi[nn];          // get neighbor for phi - 10
        //........................................................................
        nn = ijk - strideZ - 1; // neighbor index (get convention)
        m11 = Psi[nn];          // get neighbor for phi - 11
        //........................................................................
        nn = ijk + strideZ + 1; // neighbor index (get convention)
        m12 = Psi[nn];          // get neighbor for phi - 12
        //........................................................................
        nn = ijk + strideZ - 1; // neighbor index (get convention)
        m13 = Psi[nn];          // get neighbor for phi - 13
        //........................................................................
        nn = ijk - strideZ + 1; // neighbor index (get convention)
        m14 = Psi[nn];          // get neighbor for phi - 14
        //........................................................................
        nn = ijk - strideZ - strideY; // neighbor index (get convention)
        m15 = Psi[nn];                // get neighbor for phi - 15
        //........................................................................
        nn = ijk + strideZ + strideY; // neighbor index (get convention)
        m16 = Psi[nn];                // get neighbor for phi - 16
        //........................................................................
        nn = ijk + strideZ - strideY; // neighbor index (get convention)
        m17 = Psi[nn];                // get neighbor for phi - 17
        //........................................................................
        nn = ijk - strideZ + strideY; // neighbor index (get convention)
        m18 = Psi[nn];                // get neighbor for phi - 18

        psi_Laplacian =
            2.0 * 3.0 / 18.0 *
            (m1 + m2 + m3 + m4 + m5 + m6 - 6 * psi +
             0.5 * (m7 + m8 + m9 + m10 + m11 + m12 + m13 + m14 + m15 + m16 +
                    m17 + m18 - 12 * psi)); //Laplacian of electric potential
        residual_error = psi_Laplacian + rho_e / epsilon_LB;
        ResidualError[n] = residual_error;
    }
}
//extern "C" void ScaLBL_D3Q7_Poisson_ElectricField(int *neighborList, int *Map, signed char *ID, double *Psi, double *ElectricField, int SolidBC,
//        int strideY, int strideZ,int start, int finish, int Np){
//
//	int n,nn;
//    int ijk;
//    int id;
//	// distributions
//	double m1,m2,m3,m4,m5,m6,m7,m8,m9;
//	double m10,m11,m12,m13,m14,m15,m16,m17,m18;
//	double nx,ny,nz;
//
//	for (n=start; n<finish; n++){
//
//		// Get the 1D index based on regular data layout
//		ijk = Map[n];
//		//					COMPUTE THE COLOR GRADIENT
//		//........................................................................
//		//.................Read Phase Indicator Values............................
//		//........................................................................
//		nn = ijk-1;							// neighbor index (get convention)
//        id = ID[nn];
//		m1 = SolidBC==1 ? Psi[nn] : Psi[nn]*(id>0)+Psi[ijk]*(id<=0);// get neighbor for phi - 1
//		//........................................................................
//		nn = ijk+1;							// neighbor index (get convention)
//        id = ID[nn];
//		m2 = SolidBC==1 ? Psi[nn] : Psi[nn]*(id>0)+Psi[ijk]*(id<=0);// get neighbor for phi - 2
//		//........................................................................
//		nn = ijk-strideY;							// neighbor index (get convention)
//        id = ID[nn];
//		m3 = SolidBC==1 ? Psi[nn] : Psi[nn]*(id>0)+Psi[ijk]*(id<=0);// get neighbor for phi - 3
//		//........................................................................
//		nn = ijk+strideY;							// neighbor index (get convention)
//        id = ID[nn];
//		m4 = SolidBC==1 ? Psi[nn] : Psi[nn]*(id>0)+Psi[ijk]*(id<=0);// get neighbor for phi - 4
//		//........................................................................
//		nn = ijk-strideZ;						// neighbor index (get convention)
//        id = ID[nn];
//		m5 = SolidBC==1 ? Psi[nn] : Psi[nn]*(id>0)+Psi[ijk]*(id<=0);// get neighbor for phi - 5
//		//........................................................................
//		nn = ijk+strideZ;						// neighbor index (get convention)
//        id = ID[nn];
//		m6 = SolidBC==1 ? Psi[nn] : Psi[nn]*(id>0)+Psi[ijk]*(id<=0);// get neighbor for phi - 6
//		//........................................................................
//		nn = ijk-strideY-1;						// neighbor index (get convention)
//        id = ID[nn];
//		m7 = SolidBC==1 ? Psi[nn] : Psi[nn]*(id>0)+Psi[ijk]*(id<=0);// get neighbor for phi - 7
//		//........................................................................
//		nn = ijk+strideY+1;						// neighbor index (get convention)
//        id = ID[nn];
//		m8 = SolidBC==1 ? Psi[nn] : Psi[nn]*(id>0)+Psi[ijk]*(id<=0);// get neighbor for phi - 8
//		//........................................................................
//		nn = ijk+strideY-1;						// neighbor index (get convention)
//        id = ID[nn];
//		m9 = SolidBC==1 ? Psi[nn] : Psi[nn]*(id>0)+Psi[ijk]*(id<=0);// get neighbor for phi - 9
//		//........................................................................
//		nn = ijk-strideY+1;						// neighbor index (get convention)
//        id = ID[nn];
//		m10 = SolidBC==1 ? Psi[nn] : Psi[nn]*(id>0)+Psi[ijk]*(id<=0);// get neighbor for phi - 10
//		//........................................................................
//		nn = ijk-strideZ-1;						// neighbor index (get convention)
//        id = ID[nn];
//		m11 = SolidBC==1 ? Psi[nn] : Psi[nn]*(id>0)+Psi[ijk]*(id<=0);// get neighbor for phi - 11
//		//........................................................................
//		nn = ijk+strideZ+1;						// neighbor index (get convention)
//        id = ID[nn];
//		m12 = SolidBC==1 ? Psi[nn] : Psi[nn]*(id>0)+Psi[ijk]*(id<=0);// get neighbor for phi - 12
//		//........................................................................
//		nn = ijk+strideZ-1;						// neighbor index (get convention)
//        id = ID[nn];
//		m13 = SolidBC==1 ? Psi[nn] : Psi[nn]*(id>0)+Psi[ijk]*(id<=0);// get neighbor for phi - 13
//		//........................................................................
//		nn = ijk-strideZ+1;						// neighbor index (get convention)
//        id = ID[nn];
//		m14 = SolidBC==1 ? Psi[nn] : Psi[nn]*(id>0)+Psi[ijk]*(id<=0);// get neighbor for phi - 14
//		//........................................................................
//		nn = ijk-strideZ-strideY;					// neighbor index (get convention)
//        id = ID[nn];
//		m15 = SolidBC==1 ? Psi[nn] : Psi[nn]*(id>0)+Psi[ijk]*(id<=0);// get neighbor for phi - 15
//		//........................................................................
//		nn = ijk+strideZ+strideY;					// neighbor index (get convention)
//        id = ID[nn];
//		m16 = SolidBC==1 ? Psi[nn] : Psi[nn]*(id>0)+Psi[ijk]*(id<=0);// get neighbor for phi - 16
//		//........................................................................
//		nn = ijk+strideZ-strideY;					// neighbor index (get convention)
//        id = ID[nn];
//		m17 = SolidBC==1 ? Psi[nn] : Psi[nn]*(id>0)+Psi[ijk]*(id<=0);// get neighbor for phi - 17
//		//........................................................................
//		nn = ijk-strideZ+strideY;					// neighbor index (get convention)
//        id = ID[nn];
//		m18 = SolidBC==1 ? Psi[nn] : Psi[nn]*(id>0)+Psi[ijk]*(id<=0);// get neighbor for phi - 18
//		//............Compute the Color Gradient...................................
//		//nx = 1.f/6.f*(m1-m2+0.5*(m7-m8+m9-m10+m11-m12+m13-m14));
//		//ny = 1.f/6.f*(m3-m4+0.5*(m7-m8-m9+m10+m15-m16+m17-m18));
//		//nz = 1.f/6.f*(m5-m6+0.5*(m11-m12-m13+m14+m15-m16-m17+m18));
//		nx = 1.f/6.f*(m1-m2);//but looks like it needs to multiply another factor of 3
//		ny = 1.f/6.f*(m3-m4);
//		nz = 1.f/6.f*(m5-m6);
//
//		ElectricField[n] = nx;
//		ElectricField[Np+n] = ny;
//		ElectricField[2*Np+n] = nz;
//	}
//}

//extern "C" void ScaLBL_D3Q7_Poisson_getElectricField(double *dist, double *ElectricField, double tau, int Np){
//	int n;
//	// distributions
//	double f1,f2,f3,f4,f5,f6;
//	double Ex,Ey,Ez;
//    double rlx=1.0/tau;
//
//	for (n=0; n<Np; n++){
//		//........................................................................
//		// Registers to store the distributions
//		//........................................................................
//		f1 = dist[Np+n];
//		f2 = dist[2*Np+n];
//		f3 = dist[3*Np+n];
//		f4 = dist[4*Np+n];
//		f5 = dist[5*Np+n];
//		f6 = dist[6*Np+n];
//		//.................Compute the Electric Field...................................
//		//Ex = (f1-f2)*rlx*4.5;//NOTE the unit of electric field here is V/lu
//		//Ey = (f3-f4)*rlx*4.5;
//		//Ez = (f5-f6)*rlx*4.5;
//		Ex = (f1-f2)*rlx*4.0;//NOTE the unit of electric field here is V/lu
//		Ey = (f3-f4)*rlx*4.0;
//		Ez = (f5-f6)*rlx*4.0;
//		//..................Write the Electric Field.....................................
//		ElectricField[0*Np+n] = Ex;
//		ElectricField[1*Np+n] = Ey;
//		ElectricField[2*Np+n] = Ez;
//		//........................................................................
//	}
//}
