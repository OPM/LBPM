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
		TIMELOG = fopen("subphase.csv","a+");
		if (fseek(TIMELOG,0,SEEK_SET) == fseek(TIMELOG,0,SEEK_CUR))
		{
			// If timelog is empty, write a short header to list the averages
			//fprintf(TIMELOG,"--------------------------------------------------------------------------------------\n");
			fprintf(TIMELOG,"time rn rw nun nuw Fx Fy Fz iftwn ");				
			fprintf(TIMELOG,"pnc pnd pni pwc pwd pwi ");				// pressures 
			fprintf(TIMELOG,"Mwc Mwd Mwi Mnc Mnd Mni ");				// mass 
			fprintf(TIMELOG,"Pwc_x Pwd_x Pwi_x Pnc_x Pnd_x Pni_x ");	// momentum 
			fprintf(TIMELOG,"Pwc_y Pwd_y Pwi_y Pnc_y Pnd_y Pni_y ");			
			fprintf(TIMELOG,"Pwc_z Pwd_z Pwi_z Pnc_z Pnd_z Pni_z ");			
			fprintf(TIMELOG,"Kwc Kwd Kwi Knc Knd Kni ");				// kinetic energy
			fprintf(TIMELOG,"Vwc Awc Hwc Xwc ");					 	// wc region 
			fprintf(TIMELOG,"Vwd Awd Hwd Xwd ");					 	// wd region
			fprintf(TIMELOG,"Vnc Anc Hnc Xnc ");					 	// nc region
			fprintf(TIMELOG,"Vnd And Hnd Xnd ");					 	// nd region
			fprintf(TIMELOG,"Vi Ai Hi Xi\n");					 		// interface region 
			// stress tensor
		}

	}
	else{
		char LocalRankString[8];
		sprintf(LocalRankString,"%05d",Dm->rank());
		char LocalRankFilename[40];
		sprintf(LocalRankFilename,"%s%s","subphase.csv.",LocalRankString);
		TIMELOG = fopen(LocalRankFilename,"a+");
		//fprintf(TIMELOG,"--------------------------------------------------------------------------------------\n");
		fprintf(TIMELOG,"time rn rw nun nuw Fx Fy Fz iftwn ");				
		fprintf(TIMELOG,"pnc pnd pni pwc pwd pwi ");				// pressures 
		fprintf(TIMELOG,"Mwc Mwd Mwi Mnc Mnd Mni ");				// mass 
		fprintf(TIMELOG,"Pwc_x Pwd_x Pwi_x Pnc_x Pnd_x Pni_x ");	// momentum 
		fprintf(TIMELOG,"Pwc_y Pwd_y Pwi_y Pnc_y Pnd_y Pni_y ");			
		fprintf(TIMELOG,"Pwc_z Pwd_z Pwi_z Pnc_z Pnd_z Pni_z ");			
		fprintf(TIMELOG,"Kwc Kwd Kwi Knc Knd Kni ");				// kinetic energy
		fprintf(TIMELOG,"Vwc Awc Hwc Xwc ");					 	// wc region 
		fprintf(TIMELOG,"Vwd Awd Hwd Xwd ");					 	// wd region
		fprintf(TIMELOG,"Vnc Anc Hnc Xnc ");					 	// nc region
		fprintf(TIMELOG,"Vnd And Hnd Xnd ");					 	// nd region
		fprintf(TIMELOG,"Vi Ai Hi Xi\n");					 		// interface region 
	}
}


// Destructor
SubPhase::~SubPhase()
{
    if ( TIMELOG!=NULL ) { fclose(TIMELOG); }

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

void SubPhase::BulkAverage(){
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
		
	nb.reset(); wb.reset();
/*
 	//Dm->CommunicateMeshHalo(Phi);
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
	*/
	
	double nA,nB;
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
						nb.M += rho_n;						
						// velocity
						nb.Px += rho_n*nA*Vel_x(n);
						nb.Py += rho_n*nA*Vel_y(n);
						nb.Pz += rho_n*nA*Vel_z(n);

						/* // volume the excludes the interfacial region
						if (DelPhi(n) < 1e-4){
							// pressure
							pan += Pressure(n);
						}
						else{
							
						}
						*/
					}
					else{
						wb.M += rho_w;
						wb.V += 1.0;

						// velocity
						wb.Px += rho_w*nB*Vel_x(n);
						wb.Py += rho_w*nB*Vel_y(n);
						wb.Pz += rho_w*nB*Vel_z(n);

						/*
						if (DelPhi(n) < 1e-4){

							
						}
						else{
							
						}
						*/
					}
				}
			}
		}
	}
	wb.V=sumReduce( Dm->Comm, wb.V);
	wb.M=sumReduce( Dm->Comm, wb.M);
	nb.M=sumReduce( Dm->Comm, nb.M);
	wb.Px=sumReduce( Dm->Comm, wb.Px);
	wb.Py=sumReduce( Dm->Comm, wb.Py);
	wb.Pz=sumReduce( Dm->Comm, wb.Pz);
	nb.Px=sumReduce( Dm->Comm, nb.Px);
	nb.Py=sumReduce( Dm->Comm, nb.Py);
	nb.Pz=sumReduce( Dm->Comm, nb.Pz);

	if (Dm->rank() == 0){
		double saturation=wb.V/(wb.V + nb.V);
		double fractional_flow=nb.M*sqrt(wb.Px*wb.Px+wb.Py*wb.Py+wb.Pz*wb.Pz)/(wb.M*sqrt(nb.Px*nb.Px+nb.Py*nb.Py+nb.Pz*nb.Pz));
		printf("saturation = %f, fractional flow =%f \n",saturation,fractional_flow);
	}

}

inline void InterfaceTransportMeasures( double beta, double rA, double rB, double nA, double nB, 
		double nx, double ny, double nz, double ux, double uy, double uz, interface &I){
	
	double A0,A1,A2,A3,A4,A5,A6;
	double B0,B1,B2,B3,B4,B5,B6;
	double nAB,delta;
	// Instantiate mass transport distributions
	// Stationary value - distribution 0
	nAB = 1.0/(nA+nB);
	A0 = 0.3333333333333333*nA;
	B0 = 0.3333333333333333*nB;

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

void SubPhase::FullAnalysis(){
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
		
	nd.reset();	nc.reset(); wd.reset();	wc.reset();	iwn.reset();

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
	morph_n->MeasureConnectedPathway();
	nc.V = morph_n->V(); 
	nc.A = morph_n->A(); 
	nc.H = morph_n->H(); 
	nc.X = morph_n->X(); 
	// update disconnected part
	nd.V -= nc.V;
	nd.A -= nc.A;
	nd.H -= nc.H;
	nd.X -= nc.X;
	
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
	morph_w->MeasureConnectedPathway();
	wc.V = morph_w->V(); 
	wc.A = morph_w->A(); 
	wc.H = morph_w->H(); 
	wc.X = morph_w->X(); 
	// update disconnected part
	wd.V -= wc.V;
	wd.A -= wc.A;
	wd.H -= wc.H;
	wd.X -= wc.X;
	
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
					// wetting phase
					morph_i->id(i,j,k) = 0;
				}
				else {
					// non-wetting phase
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
					double phi = (rho_n-rho_w)/(rho_n+rho_w);
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
							wd.M += nB*rho_n;						
							wd.Px += nB*rho_n*ux;
							wd.Py += nB*rho_n*uy;
							wd.Pz += nB*rho_n*uz;
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
	iwn.V=sumReduce( Dm->Comm, iwn.V);
	wb.M=sumReduce( Dm->Comm, wb.M);
	nb.M=sumReduce( Dm->Comm, nb.M);
	wb.Px=sumReduce( Dm->Comm, wb.Px);
	wb.Py=sumReduce( Dm->Comm, wb.Py);
	wb.Pz=sumReduce( Dm->Comm, wb.Pz);
	nb.Px=sumReduce( Dm->Comm, nb.Px);
	nb.Py=sumReduce( Dm->Comm, nb.Py);
	nb.Pz=sumReduce( Dm->Comm, nb.Pz);

	if (Dm->rank() == 0){
		double saturation=wb.V/(wb.V + nb.V);
		double fractional_flow=nb.M*sqrt(wb.Px*wb.Px+wb.Py*wb.Py+wb.Pz*wb.Pz)/(wb.M*sqrt(nb.Px*nb.Px+nb.Py*nb.Py+nb.Pz*nb.Pz));
		printf("saturation = %f, fractional flow =%f \n",saturation,fractional_flow);
	}
}



