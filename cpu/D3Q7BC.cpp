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
		dist[iq] = -1.0*value_q + value_b*2.0/9.0;//NOTE 2/9 is the speed of sound for D3Q7 lattice
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

