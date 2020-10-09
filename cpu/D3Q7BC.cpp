// CPU Functions for D3Q7 Lattice Boltzmann Methods
// Boundary Conditions 

extern "C" void ScaLBL_Solid_Dirichlet_D3Q7(double *dist,double *BoundaryValue,int *BounceBackDist_list,int *BounceBackSolid_list,int N){

    int idx;
    int iq,ib;
    double value_b,value_q;
	for (idx=0; idx<N; idx++){
		iq = BounceBackDist_list[idx];
        ib = BounceBackSolid_list[idx];
		value_b = BoundaryValue[ib];//get boundary value from a solid site
        value_q = dist[iq];
		dist[iq] = -1.0*value_q + value_b*0.25;//NOTE 0.25 is the speed of sound for D3Q7 lattice
	}
}


extern "C" void ScaLBL_Solid_Neumann_D3Q7(double *dist,double *BoundaryValue,int *BounceBackDist_list,int *BounceBackSolid_list,int N){

    int idx;
    int iq,ib;
    double value_b,value_q;
	for (idx=0; idx<N; idx++){
		iq = BounceBackDist_list[idx];
        ib = BounceBackSolid_list[idx];
		value_b = BoundaryValue[ib];//get boundary value from a solid site
        value_q = dist[iq];
		dist[iq] = value_q + value_b;
	}
}

extern "C" void ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_z(int *list, double *dist, double Vin, int count, int Np){
	for (int idx=0; idx<count; idx++){
		int n = list[idx];
		double f0 = dist[n];
		double f1 = dist[2*Np+n];
		double f2 = dist[1*Np+n];
		double f3 = dist[4*Np+n];
		double f4 = dist[3*Np+n];
		double f6 = dist[5*Np+n];
		//...................................................
		double f5 = Vin - (f0+f1+f2+f3+f4+f6);
		dist[6*Np+n] = f5;
	}
}

extern "C" void ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_Z(int *list, double *dist, double Vout, int count, int Np){
	for (int idx=0; idx<count; idx++){
		int n = list[idx];
		double f0 = dist[n];
		double f1 = dist[2*Np+n];
		double f2 = dist[1*Np+n];
		double f3 = dist[4*Np+n];
		double f4 = dist[3*Np+n];
		double f5 = dist[6*Np+n];
		//...................................................
		double f6 = Vout - (f0+f1+f2+f3+f4+f5);
		dist[5*Np+n] = f6;
	}
}

extern "C" void ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_z(int *d_neighborList, int *list, double *dist, double Vin, int count, int Np){
    int nread,nr5;
	for (int idx=0; idx<count; idx++){
		int n = list[idx];
		double f0 = dist[n];

		nread = d_neighborList[n];
		double f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		double f3 = dist[nread];

		nread = d_neighborList[n+Np];
		double f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		double f4 = dist[nread];

		nread = d_neighborList[n+5*Np];
		double f6 = dist[nread];

		// Unknown distributions
		nr5 = d_neighborList[n+4*Np];
		double f5 = Vin - (f0+f1+f2+f3+f4+f6);
		dist[nr5] = f5;
	}
}

extern "C" void ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_Z(int *d_neighborList, int *list, double *dist, double Vout, int count, int Np){
    int nread,nr6;
	for (int idx=0; idx<count; idx++){
		int n = list[idx];
		double f0 = dist[n];

		nread = d_neighborList[n];
		double f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		double f3 = dist[nread];

		nread = d_neighborList[n+4*Np];
		double f5 = dist[nread];

		nread = d_neighborList[n+Np];
		double f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		double f4 = dist[nread];

		// unknown distributions
		nr6 = d_neighborList[n+5*Np];
		double f6 = Vout - (f0+f1+f2+f3+f4+f5);
		dist[nr6] = f6;
	}
}

extern "C" void ScaLBL_Poisson_D3Q7_BC_z(int *list, int *Map, double *Psi, double Vin, int count)
{
	int idx,n,nm;

	for (idx=0; idx<count; idx++){
		n = list[idx];
		nm = Map[n];
		Psi[nm] = Vin;
	}
}

extern "C" void ScaLBL_Poisson_D3Q7_BC_Z(int *list, int *Map, double *Psi, double Vout, int count)
{
	int idx,n,nm;

	for (idx=0; idx<count; idx++){
		n = list[idx];
		nm = Map[n];
		Psi[nm] = Vout;
	}
}

extern "C" void ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_z(int *list, double *dist, double Cin, int count, int Np){
	for (int idx=0; idx<count; idx++){
		int n = list[idx];
		double f0 = dist[n];
		double f1 = dist[2*Np+n];
		double f2 = dist[1*Np+n];
		double f3 = dist[4*Np+n];
		double f4 = dist[3*Np+n];
		double f6 = dist[5*Np+n];
		//...................................................
		double f5 = Cin - (f0+f1+f2+f3+f4+f6);
		dist[6*Np+n] = f5;
	}
}

extern "C" void ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_Z(int *list, double *dist, double Cout, int count, int Np){
	for (int idx=0; idx<count; idx++){
		int n = list[idx];
		double f0 = dist[n];
		double f1 = dist[2*Np+n];
		double f2 = dist[1*Np+n];
		double f3 = dist[4*Np+n];
		double f4 = dist[3*Np+n];
		double f5 = dist[6*Np+n];
		//...................................................
		double f6 = Cout - (f0+f1+f2+f3+f4+f5);
		dist[5*Np+n] = f6;
	}
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_z(int *d_neighborList, int *list, double *dist, double Cin, int count, int Np){
    int nread,nr5;
	for (int idx=0; idx<count; idx++){
		int n = list[idx];
		double f0 = dist[n];

		nread = d_neighborList[n];
		double f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		double f3 = dist[nread];

		nread = d_neighborList[n+Np];
		double f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		double f4 = dist[nread];

		nread = d_neighborList[n+5*Np];
		double f6 = dist[nread];

		// Unknown distributions
		nr5 = d_neighborList[n+4*Np];
		double f5 = Cin - (f0+f1+f2+f3+f4+f6);
		dist[nr5] = f5;
	}
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_Z(int *d_neighborList, int *list, double *dist, double Cout, int count, int Np){
    int nread,nr6;
	for (int idx=0; idx<count; idx++){
		int n = list[idx];
		double f0 = dist[n];

		nread = d_neighborList[n];
		double f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		double f3 = dist[nread];

		nread = d_neighborList[n+4*Np];
		double f5 = dist[nread];

		nread = d_neighborList[n+Np];
		double f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		double f4 = dist[nread];

		// unknown distributions
		nr6 = d_neighborList[n+5*Np];
		double f6 = Cout - (f0+f1+f2+f3+f4+f5);
		dist[nr6] = f6;
	}
}
