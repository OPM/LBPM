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
	
	for (k=kmin; k<kmax; k++){
		for (j=jmin; j<Ny-1; j++){
			for (i=imin; i<Nx-1; i++){
				n = k*Nx*Ny + j*Nx + i;
				// Compute volume averages
				if ( Dm->id[n] > 0 ){
					// compute density
					double rho_n = Rho_n(n);
					double rho_w = Rho_w(n);
					double phi = (rho_n-rho_w)/(rho_n+rho_w);
					Phi(n) = phi;
					
					if ( phi > 0.0 ){
						nb.V += 1.0;
						nb.M += rho_n;						
						// velocity
						nb.Px += rho_n*Vel_x(n);
						nb.Py += rho_n*Vel_y(n);
						nb.Pz += rho_n*Vel_z(n);

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
						wb.Px += rho_w*Vel_x(n);
						wb.Py += rho_w*Vel_y(n);
						wb.Pz += rho_w*Vel_z(n);

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
