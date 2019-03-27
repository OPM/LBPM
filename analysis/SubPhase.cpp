#include "analysis/SubPhase.h"

// Constructor
SubPhase::SubPhase(std::shared_ptr <Domain> dm):
	Dm(dm)
{
	Nx=dm->Nx; Ny=dm->Ny; Nz=dm->Nz;
	Volume=(Nx-2)*(Ny-2)*(Nz-2)*Dm->nprocx()*Dm->nprocy()*Dm->nprocz()*1.0;
	
	morph_w = std::shared_ptr<Minkowski>(new Minkowski(Dm));
	morph_n = std::shared_ptr<Minkowski>(new Minkowski(Dm));
	morph_i = std::shared_ptr<Minkowski>(new Minkowski(Dm));

	// Global arrays
	PhaseID.resize(Nx,Ny,Nz);       PhaseID.fill(0);
	Label_WP.resize(Nx,Ny,Nz);      Label_WP.fill(0);
	Label_NWP.resize(Nx,Ny,Nz);     Label_NWP.fill(0);
	Rho_n.resize(Nx,Ny,Nz);       	Rho_n.fill(0);
	Rho_w.resize(Nx,Ny,Nz);       	Rho_w.fill(0);
	Pressure.resize(Nx,Ny,Nz);      Pressure.fill(0);
	Phi.resize(Nx,Ny,Nz);         	Phi.fill(0);
	DelPhi.resize(Nx,Ny,Nz);        DelPhi.fill(0);
	Vel_x.resize(Nx,Ny,Nz);         Vel_x.fill(0);	    // Gradient of the phase indicator field
	Vel_y.resize(Nx,Ny,Nz);         Vel_y.fill(0);
	Vel_z.resize(Nx,Ny,Nz);         Vel_z.fill(0);
	SDs.resize(Nx,Ny,Nz);         	SDs.fill(0);
	//.........................................

	//.........................................
	if (Dm->rank()==0){
		bool WriteHeader=false;
		SUBPHASE = fopen("subphase.csv","r");
		if (SUBPHASE != NULL)
			fclose(SUBPHASE);
		else
			WriteHeader=true;

		SUBPHASE = fopen("subphase.csv","a+");
		if (WriteHeader)
		{
			// If timelog is empty, write a short header to list the averages
			//fprintf(SUBPHASE,"--------------------------------------------------------------------------------------\n");
			fprintf(SUBPHASE,"time rn rw nun nuw Fx Fy Fz iftwn ");				
			fprintf(SUBPHASE,"pwc pwd pnc pnd ");						// pressures 
			fprintf(SUBPHASE,"Mwc Mwd Mwi Mnc Mnd Mni ");				// mass 
			fprintf(SUBPHASE,"Pwc_x Pwd_x Pwi_x Pnc_x Pnd_x Pni_x ");	// momentum 
			fprintf(SUBPHASE,"Pwc_y Pwd_y Pwi_y Pnc_y Pnd_y Pni_y ");			
			fprintf(SUBPHASE,"Pwc_z Pwd_z Pwi_z Pnc_z Pnd_z Pni_z ");			
			fprintf(SUBPHASE,"Kwc Kwd Kwi Knc Knd Kni ");				// kinetic energy
			fprintf(SUBPHASE,"Vwc Awc Hwc Xwc ");					 	// wc region 
			fprintf(SUBPHASE,"Vwd Awd Hwd Xwd Nwd ");					// wd region
			fprintf(SUBPHASE,"Vnc Anc Hnc Xnc ");					 	// nc region
			fprintf(SUBPHASE,"Vnd And Hnd Xnd Nnd ");					// nd region
			fprintf(SUBPHASE,"Vi Ai Hi Xi ");					 		// interface region 
			fprintf(SUBPHASE,"Vic Aic Hic Xic Nic\n");					// interface region 

			// stress tensor?
		}

	}
	else{
		char LocalRankString[8];
		sprintf(LocalRankString,"%05d",Dm->rank());
		char LocalRankFilename[40];
		sprintf(LocalRankFilename,"%s%s","subphase.csv.",LocalRankString);
		SUBPHASE = fopen(LocalRankFilename,"a+");
		//fprintf(SUBPHASE,"--------------------------------------------------------------------------------------\n");
		fprintf(SUBPHASE,"time rn rw nun nuw Fx Fy Fz iftwn ");				
		fprintf(SUBPHASE,"pwc pwd pnc pnd ");						// pressures 
		fprintf(SUBPHASE,"Mwc Mwd Mwi Mnc Mnd Mni ");				// mass 
		fprintf(SUBPHASE,"Pwc_x Pwd_x Pwi_x Pnc_x Pnd_x Pni_x ");	// momentum 
		fprintf(SUBPHASE,"Pwc_y Pwd_y Pwi_y Pnc_y Pnd_y Pni_y ");			
		fprintf(SUBPHASE,"Pwc_z Pwd_z Pwi_z Pnc_z Pnd_z Pni_z ");			
		fprintf(SUBPHASE,"Kwc Kwd Kwi Knc Knd Kni ");				// kinetic energy
		fprintf(SUBPHASE,"Vwc Awc Hwc Xwc ");					 	// wc region 
		fprintf(SUBPHASE,"Vwd Awd Hwd Xwd Nwd ");					// wd region
		fprintf(SUBPHASE,"Vnc Anc Hnc Xnc ");					 	// nc region
		fprintf(SUBPHASE,"Vnd And Hnd Xnd Nnd ");					// nd region
		fprintf(SUBPHASE,"Vi Ai Hi Xi ");					 		// interface region 
		fprintf(SUBPHASE,"Vic Aic Hic Xic Nic\n");					// interface region 
	}
	
	if (Dm->rank()==0){
		bool WriteHeader=false;
		TIMELOG = fopen("timelog.csv","r");
		if (TIMELOG != NULL)
			fclose(TIMELOG);
		else
			WriteHeader=true;

		TIMELOG = fopen("timelog.csv","a+");
		if (WriteHeader)
		{
			// If timelog is empty, write a short header to list the averages
			//fprintf(TIMELOG,"--------------------------------------------------------------------------------------\n");
			fprintf(TIMELOG,"sw krw krn qw qn pw pn\n");				
		}
	}
}


// Destructor
SubPhase::~SubPhase()
{
    if ( SUBPHASE!=NULL ) { fclose(SUBPHASE); }

}

void SubPhase::Write(int timestep)
{
	if (Dm->rank()==0){
		fprintf(SUBPHASE,"%i %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g ",timestep,rho_n,rho_w,nu_n,nu_w,Fx,Fy,Fz,gamma_wn); 
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g ",gwc.p, gwd.p, gnc.p, gnd.p);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g %.5g %.5g ",gwc.M, gwd.M, giwn.Mw, gnc.M, gnd.M, giwn.Mn);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g %.5g %.5g ",gwc.Px, gwd.Px, giwn.Px, gnc.Px, gnd.Px, giwn.Px);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g %.5g %.5g ",gwc.Py, gwd.Py, giwn.Py, gnc.Py, gnd.Py, giwn.Py);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g %.5g %.5g ",gwc.Pz, gwd.Pz, giwn.Pz, gnc.Pz, gnd.Pz, giwn.Pz);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g %.5g %.5g ",gwc.K, gwd.K, giwn.Kw, gnc.K, gnd.K, giwn.Kn);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g ",gwc.V, gwc.A, gwc.H, gwc.X);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g %i ",gwd.V, gwd.A, gwd.H, gwd.X, gwd.Nc);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g ",gnc.V, gnc.A, gnc.H, gnc.X);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g %i ",gnd.V, gnd.A, gnd.H, gnd.X, gnd.Nc);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g ",giwn.V, giwn.A, giwn.H, giwn.X);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g %i\n",giwnc.V, giwnc.A, giwnc.H, giwnc.X, giwnc.Nc);
		fflush(SUBPHASE);
	}
	else{
		fprintf(SUBPHASE,"%i %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g ",timestep,rho_n,rho_w,nu_n,nu_w,Fx,Fy,Fz,gamma_wn);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g ",wc.p, wd.p, nc.p, nd.p);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g %.5g %.5g ",wc.M, wd.M, iwn.Mw, nc.M, nd.M, iwn.Mn);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g %.5g %.5g ",wc.Px, wd.Px, iwn.Px, nc.Px, nd.Px, iwn.Px);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g %.5g %.5g ",wc.Py, wd.Py, iwn.Py, nc.Py, nd.Py, iwn.Py);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g %.5g %.5g ",wc.Pz, wd.Pz, iwn.Pz, nc.Pz, nd.Pz, iwn.Pz);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g %.5g %.5g ",wc.K, wd.K, iwn.Kw, nc.K, nd.K, iwn.Kn);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g ",wc.V, wc.A, wc.H, wc.X);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g %i ",wd.V, wd.A, wd.H, wd.X, wd.Nc);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g ",nc.V, nc.A, nc.H, nc.X);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g %i ",nd.V, nd.A, nd.H, nd.X, nd.Nc);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g ",iwn.V, iwn.A, iwn.H, iwn.X);
		fprintf(SUBPHASE,"%.5g %.5g %.5g %.5g\n",iwnc.V, iwnc.A, iwnc.H, iwnc.X);
	}

}

void SubPhase::SetParams(double rhoA, double rhoB, double tauA, double tauB, double force_x, double force_y, double force_z, double alpha, double B)
{
	Fx = force_x;
	Fy = force_y;
	Fz = force_z;
	rho_n = rhoA;
	rho_w = rhoB;
	nu_n = (tauA-0.5)/3.f;
	nu_w = (tauB-0.5)/3.f;
	gamma_wn = 5.796*alpha;
	beta = B;
}

void SubPhase::Basic(){
	int i,j,k,n,imin,jmin,kmin,kmax;

	// If external boundary conditions are set, do not average over the inlet
	kmin=1; kmax=Nz-1;
	if (Dm->BoundaryCondition > 0 && Dm->kproc() == 0) kmin=4;
	if (Dm->BoundaryCondition > 0 && Dm->kproc() == Dm->nprocz()-1) kmax=Nz-4;

	imin=jmin=1;
	// If inlet/outlet layers exist use these as default
	if (Dm->inlet_layers_x > 0) imin = Dm->inlet_layers_x;
	if (Dm->inlet_layers_y > 0) jmin = Dm->inlet_layers_y;
	if (Dm->inlet_layers_z > 0) kmin = Dm->inlet_layers_z;
	if (Dm->outlet_layers_z > 0) kmax = Dm->outlet_layers_z;
	
	nb.reset(); wb.reset();

	double nA,nB;
	double count_w = 0.0;
	double count_n = 0.0;
	for (k=kmin; k<kmax; k++){
		for (j=jmin; j<Ny-1; j++){
			for (i=imin; i<Nx-1; i++){
				n = k*Nx*Ny + j*Nx + i;
				// Compute volume averages
				if ( Dm->id[n] > 0 ){
					// compute density
					double nA = Rho_n(n);
					double nB = Rho_w(n);
					double phi = (nA-nB)/(nA+nB);
					Phi(n) = phi;
					
					if ( phi > 0.0 ){
						nb.V += 1.0;
						nb.M += nA*rho_n;						
						// velocity
						nb.Px += rho_n*nA*Vel_x(n);
						nb.Py += rho_n*nA*Vel_y(n);
						nb.Pz += rho_n*nA*Vel_z(n);
					}
					else{
						wb.M += nB*rho_w;
						wb.V += 1.0;

						// velocity
						wb.Px += rho_w*nB*Vel_x(n);
						wb.Py += rho_w*nB*Vel_y(n);
						wb.Pz += rho_w*nB*Vel_z(n);
					}
					if ( phi > 0.99 ){
						nb.p += Pressure(n);
						count_n += 1.0;
					}
					else if ( phi < -0.99 ){
						wb.p += Pressure(n);
						count_w += 1.0;
					}
				}
			}
		}
	}
	gwb.V=sumReduce( Dm->Comm, wb.V);
	gnb.V=sumReduce( Dm->Comm, nb.V);
	gwb.M=sumReduce( Dm->Comm, wb.M);
	gnb.M=sumReduce( Dm->Comm, nb.M);
	gwb.Px=sumReduce( Dm->Comm, wb.Px);
	gwb.Py=sumReduce( Dm->Comm, wb.Py);
	gwb.Pz=sumReduce( Dm->Comm, wb.Pz);
	gnb.Px=sumReduce( Dm->Comm, nb.Px);
	gnb.Py=sumReduce( Dm->Comm, nb.Py);
	gnb.Pz=sumReduce( Dm->Comm, nb.Pz);
	
	count_w=sumReduce( Dm->Comm, count_w);
	count_n=sumReduce( Dm->Comm, count_n);
	gwb.p=sumReduce( Dm->Comm, wb.p) / count_w;
	gnb.p=sumReduce( Dm->Comm, nb.p) / count_n;

	if (Dm->rank() == 0){
		double force_mag = sqrt(Fx*Fx+Fy*Fy+Fz*Fz);
		double dir_x = Fx/force_mag;
		double dir_y = Fy/force_mag;
		double dir_z = Fz/force_mag;
		if (force_mag == 0.0){
			// default to z direction
			dir_z = 1.0;
			force_mag = 1.0;
		}
		double saturation=gwb.V/(gwb.V + gnb.V);
		double water_flow_rate=gwb.V*(gwb.Px*dir_x + gwb.Py*dir_y + gwb.Pz*dir_z)/gwb.M;
		double not_water_flow_rate=gnb.V*(gnb.Px*dir_x + gnb.Py*dir_y + gnb.Pz*dir_z)/gnb.M;
		double total_flow_rate = water_flow_rate + not_water_flow_rate;
		double fractional_flow= water_flow_rate / total_flow_rate;
		
		double krn = nu_n*not_water_flow_rate / force_mag;
		double krw = nu_w*water_flow_rate / force_mag;
		//printf("   water saturation = %f, fractional flow =%f \n",saturation,fractional_flow);
		fprintf(TIMELOG,"%.5g %.5g %.5g %.5g %.5g %.5g %.5g\n",saturation,krw,krn,water_flow_rate,not_water_flow_rate, gwb.p, gnb.p); 
		fflush(TIMELOG);
	}

}

inline void InterfaceTransportMeasures( double beta, double rA, double rB, double nA, double nB, 
		double nx, double ny, double nz, double ux, double uy, double uz, interface &I){
	
	double A1,A2,A3,A4,A5,A6;
	double B1,B2,B3,B4,B5,B6;
	double nAB,delta;
	// Instantiate mass transport distributions
	// Stationary value - distribution 0
	nAB = 1.0/(nA+nB);
	//...............................................
	// q = 0,2,4
	// Cq = {1,0,0}, {0,1,0}, {0,0,1}
	delta = beta*nA*nB*nAB*0.1111111111111111*nx;
	if (!(nA*nB*nAB>0)) delta=0;
	A1 = nA*(0.1111111111111111*(1+4.5*ux))+delta;
	B1 = nB*(0.1111111111111111*(1+4.5*ux))-delta;
	A2 = nA*(0.1111111111111111*(1-4.5*ux))-delta;
	B2 = nB*(0.1111111111111111*(1-4.5*ux))+delta;

	//...............................................
	// Cq = {0,1,0}
	delta = beta*nA*nB*nAB*0.1111111111111111*ny;
	if (!(nA*nB*nAB>0)) delta=0;
	A3 = nA*(0.1111111111111111*(1+4.5*uy))+delta;
	B3 = nB*(0.1111111111111111*(1+4.5*uy))-delta;
	A4 = nA*(0.1111111111111111*(1-4.5*uy))-delta;
	B4 = nB*(0.1111111111111111*(1-4.5*uy))+delta;

	//...............................................
	// q = 4
	// Cq = {0,0,1}
	delta = beta*nA*nB*nAB*0.1111111111111111*nz;
	if (!(nA*nB*nAB>0)) delta=0;
	A5 = nA*(0.1111111111111111*(1+4.5*uz))+delta;
	B5 = nB*(0.1111111111111111*(1+4.5*uz))-delta;
	A6 = nA*(0.1111111111111111*(1-4.5*uz))-delta;
	B6 = nB*(0.1111111111111111*(1-4.5*uz))+delta;
	
	double unx = (A1-A2);
	double uny = (A3-A4);
	double unz = (A5-A6);
	double uwx = (B1-B2);
	double uwy = (B3-B4);
	double uwz = (B5-B6);
	
	I.Mn += rA*nA;
	I.Mw += rB*nB;
	I.Pnx += rA*nA*unx;
	I.Pny += rA*nA*uny;
	I.Pnz += rA*nA*unz;
	I.Pwx += rB*nB*uwx;
	I.Pwy += rB*nB*uwy;
	I.Pwz += rB*nB*uwz;
	I.Kn += rA*nA*(unx*unx + uny*uny + unz*unz);
	I.Kw += rB*nB*(uwx*uwx + uwy*uwy + uwz*uwz);

}

void SubPhase::Full(){
	int i,j,k,n,imin,jmin,kmin,kmax;

	// If external boundary conditions are set, do not average over the inlet
	kmin=1; kmax=Nz-1;
	if (Dm->BoundaryCondition > 0 && Dm->kproc() == 0) kmin=4;
	if (Dm->BoundaryCondition > 0 && Dm->kproc() == Dm->nprocz()-1) kmax=Nz-4;

	imin=jmin=1;
	// If inlet layers exist use these as default
	if (Dm->inlet_layers_x > 0) imin = Dm->inlet_layers_x;
	if (Dm->inlet_layers_y > 0) jmin = Dm->inlet_layers_y;
	if (Dm->inlet_layers_z > 0) kmin = Dm->inlet_layers_z;
		
	nd.reset();	nc.reset(); wd.reset();	wc.reset();	iwn.reset();	iwnc.reset();

 	Dm->CommunicateMeshHalo(Phi);
	for (int k=1; k<Nz-1; k++){
		for (int j=1; j<Ny-1; j++){
			for (int i=1; i<Nx-1; i++){
				// Compute all of the derivatives using finite differences
				double fx = 0.5*(Phi(i+1,j,k) - Phi(i-1,j,k));
				double fy = 0.5*(Phi(i,j+1,k) - Phi(i,j-1,k));
				double fz = 0.5*(Phi(i,j,k+1) - Phi(i,j,k-1));
				DelPhi(i,j,k) = sqrt(fx*fx+fy*fy+fz*fz);
			}
		}
	}
 	Dm->CommunicateMeshHalo(DelPhi);

 	/*  Set up geometric analysis of each region */
	
 	// non-wetting
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				n = k*Nx*Ny+j*Nx+i;
				if (!(Dm->id[n] > 0)){
					// Solid phase
					morph_n->id(i,j,k) = 1;

				}
				else if (Phi(n) > 0.0){
					// non-wetting phase
					morph_n->id(i,j,k) = 0;
				}
				else {
					// wetting phase
					morph_n->id(i,j,k) = 1;
				}
			}
		}
	}
	// measure the whole object
	morph_n->MeasureObject();
	nd.V = morph_n->V(); 
	nd.A = morph_n->A(); 
	nd.H = morph_n->H(); 
	nd.X = morph_n->X(); 
	// measure only the connected part
	nd.Nc = morph_n->MeasureConnectedPathway();
	nc.V = morph_n->V(); 
	nc.A = morph_n->A(); 
	nc.H = morph_n->H(); 
	nc.X = morph_n->X(); 
	// update disconnected part
	nd.V -= nc.V;
	nd.A -= nc.A;
	nd.H -= nc.H;
	nd.X -= nc.X;

	// compute global entities
	gnc.V=sumReduce( Dm->Comm, nc.V);
	gnc.A=sumReduce( Dm->Comm, nc.A);
	gnc.H=sumReduce( Dm->Comm, nc.H);
	gnc.X=sumReduce( Dm->Comm, nc.X);
	gnd.V=sumReduce( Dm->Comm, nd.V);
	gnd.A=sumReduce( Dm->Comm, nd.A);
	gnd.H=sumReduce( Dm->Comm, nd.H);
	gnd.X=sumReduce( Dm->Comm, nd.X);
	gnd.Nc = nd.Nc;
 	// wetting
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				n = k*Nx*Ny+j*Nx+i;
				if (!(Dm->id[n] > 0)){
					// Solid phase
					morph_w->id(i,j,k) = 1;

				}
				else if (Phi(n) < 0.0){
					// wetting phase
					morph_w->id(i,j,k) = 0;
				}
				else {
					// non-wetting phase
					morph_w->id(i,j,k) = 1;
				}
			}
		}
	}	
	morph_w->MeasureObject();
	wd.V = morph_w->V(); 
	wd.A = morph_w->A(); 
	wd.H = morph_w->H(); 
	wd.X = morph_w->X(); 
	// measure only the connected part
	wd.Nc = morph_w->MeasureConnectedPathway();
	wc.V = morph_w->V(); 
	wc.A = morph_w->A(); 
	wc.H = morph_w->H(); 
	wc.X = morph_w->X(); 
	// update disconnected part
	wd.V -= wc.V;
	wd.A -= wc.A;
	wd.H -= wc.H;
	wd.X -= wc.X;
	// compute global entities
	gwc.V=sumReduce( Dm->Comm, wc.V);
	gwc.A=sumReduce( Dm->Comm, wc.A);
	gwc.H=sumReduce( Dm->Comm, wc.H);
	gwc.X=sumReduce( Dm->Comm, wc.X);
	gwd.V=sumReduce( Dm->Comm, wd.V);
	gwd.A=sumReduce( Dm->Comm, wd.A);
	gwd.H=sumReduce( Dm->Comm, wd.H);
	gwd.X=sumReduce( Dm->Comm, wd.X);
	gwd.Nc = wd.Nc;
	
 	/*  Set up geometric analysis of interface region */
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				n = k*Nx*Ny+j*Nx+i;
				if (!(Dm->id[n] > 0)){
					// Solid phase
					morph_i->id(i,j,k) = 1;
				}
				else if (DelPhi(n) > 1e-4){
					// interface
					morph_i->id(i,j,k) = 0;
				}
				else {
					// not interface
					morph_i->id(i,j,k) = 1;
				}
			}
		}
	}	
	morph_i->MeasureObject();
	iwn.V = morph_i->V(); 
	iwn.A = morph_i->A(); 
	iwn.H = morph_i->H(); 
	iwn.X = morph_i->X(); 
	giwn.V=sumReduce( Dm->Comm, iwn.V);
	giwn.A=sumReduce( Dm->Comm, iwn.A);
	giwn.H=sumReduce( Dm->Comm, iwn.H);
	giwn.X=sumReduce( Dm->Comm, iwn.X);
	// measure only the connected part
	iwnc.Nc = morph_i->MeasureConnectedPathway();
	iwnc.V = morph_i->V(); 
	iwnc.A = morph_i->A(); 
	iwnc.H = morph_i->H(); 
	iwnc.X = morph_i->X(); 
	giwnc.V=sumReduce( Dm->Comm, iwnc.V);
	giwnc.A=sumReduce( Dm->Comm, iwnc.A);
	giwnc.H=sumReduce( Dm->Comm, iwnc.H);
	giwnc.X=sumReduce( Dm->Comm, iwnc.X);
	giwnc.Nc = iwnc.Nc;

	double vol_nc_bulk = 0.0;
	double vol_wc_bulk = 0.0;
	double vol_nd_bulk = 0.0;
	double vol_wd_bulk = 0.0;
	for (k=kmin; k<kmax; k++){
		for (j=jmin; j<Ny-1; j++){
			for (i=imin; i<Nx-1; i++){
				n = k*Nx*Ny + j*Nx + i;
				// Compute volume averages
				if ( Dm->id[n] > 0 ){
					// compute density
					double nA = Rho_n(n);
					double nB = Rho_w(n);
					double phi = (nA-nB)/(nA+nB);
					double ux = Vel_x(n);
					double uy = Vel_y(n);
					double uz = Vel_z(n);
					Phi(n) = phi;

					if (DelPhi(n) > 1e-4){
						// interface region
						double nx = 0.5*(Phi(i+1,j,k)-Phi(i-1,j,k));
						double ny = 0.5*(Phi(i,j+1,k)-Phi(i,j-1,k));
						double nz = 0.5*(Phi(i,j,k+1)-Phi(i,j,k-1));
						InterfaceTransportMeasures(  beta,  rho_w,  rho_n,  nA, nB, nx, ny, nz, ux, uy, uz, iwn);
					}
					else if ( phi > 0.0){
						if (morph_n->label(i,j,k) > 0 ){
							vol_nd_bulk += 1.0;
							nd.M += nA*rho_n;						
							nd.Px += nA*rho_n*ux;
							nd.Py += nA*rho_n*uy;
							nd.Pz += nA*rho_n*uz;
							nd.K += nA*rho_n*(ux*ux + uy*uy + uz*uz);
							nd.p += Pressure(n);
						}
						else{
							vol_nc_bulk += 1.0;
							nc.M += nA*rho_n;						
							nc.Px += nA*rho_n*ux;
							nc.Py += nA*rho_n*uy;
							nc.Pz += nA*rho_n*uz;
							nc.K += nA*rho_n*(ux*ux + uy*uy + uz*uz);
							nc.p += Pressure(n);
						}
					}
					else{
						// water region
						if (morph_w->label(i,j,k) > 0 ){
							vol_wd_bulk += 1.0;
							wd.M += nB*rho_w;						
							wd.Px += nB*rho_w*ux;
							wd.Py += nB*rho_w*uy;
							wd.Pz += nB*rho_w*uz;
							wd.K += nB*rho_w*(ux*ux + uy*uy + uz*uz);
							wd.p += Pressure(n);
						}
						else{
							vol_wc_bulk += 1.0;
							wc.M += nB*rho_w;						
							wc.Px += nB*rho_w*ux;
							wc.Py += nB*rho_w*uy;
							wc.Pz += nB*rho_w*uz;
							wc.K += nB*rho_w*(ux*ux + uy*uy + uz*uz);
							wc.p += Pressure(n);
						}
					}
				}
			}
		}
	}
	gnd.M=sumReduce( Dm->Comm, nd.M);
	gnd.Px=sumReduce( Dm->Comm, nd.Px);
	gnd.Py=sumReduce( Dm->Comm, nd.Py);
	gnd.Pz=sumReduce( Dm->Comm, nd.Pz);
	gnd.K=sumReduce( Dm->Comm, nd.K);
	gnd.p=sumReduce( Dm->Comm, nd.p);

	gwd.M=sumReduce( Dm->Comm, wd.M);
	gwd.Px=sumReduce( Dm->Comm, wd.Px);
	gwd.Py=sumReduce( Dm->Comm, wd.Py);
	gwd.Pz=sumReduce( Dm->Comm, wd.Pz);
	gwd.K=sumReduce( Dm->Comm, wd.K);
	gwd.p=sumReduce( Dm->Comm, wd.p);
	
	gnc.M=sumReduce( Dm->Comm, nc.M);
	gnc.Px=sumReduce( Dm->Comm, nc.Px);
	gnc.Py=sumReduce( Dm->Comm, nc.Py);
	gnc.Pz=sumReduce( Dm->Comm, nc.Pz);
	gnc.K=sumReduce( Dm->Comm, nc.K);
	gnc.p=sumReduce( Dm->Comm, nc.p);

	gwc.M=sumReduce( Dm->Comm, wc.M);
	gwc.Px=sumReduce( Dm->Comm, wc.Px);
	gwc.Py=sumReduce( Dm->Comm, wc.Py);
	gwc.Pz=sumReduce( Dm->Comm, wc.Pz);
	gwc.K=sumReduce( Dm->Comm, wc.K);
	gwc.p=sumReduce( Dm->Comm, wc.p);
	
	giwn.Mn=sumReduce( Dm->Comm, iwn.Mn);
	giwn.Pnx=sumReduce( Dm->Comm, iwn.Pnx);
	giwn.Pny=sumReduce( Dm->Comm, iwn.Pny);
	giwn.Pnz=sumReduce( Dm->Comm, iwn.Pnz);
	giwn.Kn=sumReduce( Dm->Comm, iwn.Kn);
	giwn.Mw=sumReduce( Dm->Comm, iwn.Mw);
	giwn.Pwx=sumReduce( Dm->Comm, iwn.Pwx);
	giwn.Pwy=sumReduce( Dm->Comm, iwn.Pwy);
	giwn.Pwz=sumReduce( Dm->Comm, iwn.Pwz);
	giwn.Kw=sumReduce( Dm->Comm, iwn.Kw);
	
	// pressure averaging
	if (vol_wc_bulk > 0.0)
		wc.p = wc.p /vol_wc_bulk;
	if (vol_nc_bulk > 0.0)
		nc.p = nc.p /vol_nc_bulk;
	if (vol_wd_bulk > 0.0)
		wd.p = wd.p /vol_wd_bulk;
	if (vol_nd_bulk > 0.0)
		nd.p = nd.p /vol_nd_bulk;

	vol_wc_bulk=sumReduce( Dm->Comm, vol_wc_bulk);
	vol_wd_bulk=sumReduce( Dm->Comm, vol_wd_bulk);
	vol_nc_bulk=sumReduce( Dm->Comm, vol_nc_bulk);
	vol_nd_bulk=sumReduce( Dm->Comm, vol_nd_bulk);
	
	if (vol_wc_bulk > 0.0)
		gwc.p = gwc.p /vol_wc_bulk;
	if (vol_nc_bulk > 0.0)
		gnc.p = gnc.p /vol_nc_bulk;
	if (vol_wd_bulk > 0.0)
		gwd.p = gwd.p /vol_wd_bulk;
	if (vol_nd_bulk > 0.0)
		gnd.p = gnd.p /vol_nd_bulk;
}



