#include "analysis/SubPhase.h"

// Constructor
SubPhase::SubPhase(std::shared_ptr <Domain> dm):
	n_nw_pts(0), n_ns_pts(0), n_ws_pts(0), n_nws_pts(0), n_local_sol_pts(0), n_local_nws_pts(0),
    n_nw_tris(0), n_ns_tris(0), n_ws_tris(0), n_nws_seg(0), n_local_sol_tris(0),
    nc(0), kstart(0), kfinish(0), fluid_isovalue(0), solid_isovalue(0),	Volume(0),
    TIMELOG(NULL), NWPLOG(NULL), WPLOG(NULL), 
    Dm(dm), NumberComponents_WP(0), NumberComponents_NWP(0), trimdist(0),
    porosity(0), poreVol(0), awn(0), ans(0), aws(0), lwns(0), wp_volume(0), nwp_volume(0),
    As(0), dummy(0), vol_w(0), vol_n(0), sat_w(0), sat_w_previous(0),
    pan(0), paw(0), pan_global(0), paw_global(0), vol_w_global(0), vol_n_global(0),
    awn_global(0), ans_global(0),aws_global(0), lwns_global(0), efawns(0), efawns_global(0),
    Jwn(0), Jwn_global(0), Kwn(0), Kwn_global(0), KNwns(0), KNwns_global(0),
    KGwns(0), KGwns_global(0), trawn(0), trawn_global(0), trJwn(0), trJwn_global(0),
    trRwn(0), trRwn_global(0), nwp_volume_global(0), wp_volume_global(0),
    As_global(0), wwndnw_global(0), wwnsdnwn_global(0), Jwnwwndnw_global(0), dEs(0), dAwn(0), dAns(0)
{
	Nx=dm->Nx; Ny=dm->Ny; Nz=dm->Nz;
	Volume=(Nx-2)*(Ny-2)*(Nz-2)*Dm->nprocx()*Dm->nprocy()*Dm->nprocz()*1.0;

	TempID = new char[Nx*Ny*Nz];
	
	morph_w = std::shared_ptr<Minkowski>(new Minkowski(Dm));
	morph_n = std::shared_ptr<Minkowski>(new Minkowski(Dm));
	morph_i = std::shared_ptr<Minkowski>(new Minkowski(Dm));

	// Global arrays
	PhaseID.resize(Nx,Ny,Nz);       PhaseID.fill(0);
	Label_WP.resize(Nx,Ny,Nz);      Label_WP.fill(0);
	Label_NWP.resize(Nx,Ny,Nz);     Label_NWP.fill(0);
	Density.resize(Nx,Ny,Nz);       Density.fill(0);
	Pressure.resize(Nx,Ny,Nz);      Pressure.fill(0);
	Phi.resize(Nx,Ny,Nz);         	Phi.fill(0);
	DelPhi.resize(Nx,Ny,Nz);        DelPhi.fill(0);
	Vel_x.resize(Nx,Ny,Nz);         Vel_x.fill(0);	    // Gradient of the phase indicator field
	Vel_y.resize(Nx,Ny,Nz);         Vel_y.fill(0);
	Vel_z.resize(Nx,Ny,Nz);         Vel_z.fill(0);
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
			fprintf(TIMELOG,"Vi Ai Hi Xi ");					 		// interface region 

			// stress tensor
			
			fprintf(TIMELOG,"wwndnw wwnsdnwn Jwnwwndnw "); 	//kinematic quantities,
			fprintf(TIMELOG,"Vw Aw Jw Xw "); 			//miknowski measures,
			fprintf(TIMELOG,"Vn An Jn Xn\n"); 			//miknowski measures,
		}

	}
	else{
		char LocalRankString[8];
		sprintf(LocalRankString,"%05d",Dm->rank());
		char LocalRankFilename[40];
		sprintf(LocalRankFilename,"%s%s","subphase.csv.",LocalRankString);
		TIMELOG = fopen(LocalRankFilename,"a+");
		//fprintf(TIMELOG,"--------------------------------------------------------------------------------------\n");
		fprintf(TIMELOG,"time rn rw nun nuw Fx Fy Fz iftwn ");;								// Timestep, 
		fprintf(TIMELOG,"sw pw pn awn ans aws Jwn Kwn lwns cwns KNwns KGwns ");	// Scalar averages
		fprintf(TIMELOG,"vawx vawy vawz vanx vany vanz ");			// Velocity averages
		fprintf(TIMELOG,"vawnx vawny vawnz vawnsx vawnsy vawnsz ");
		fprintf(TIMELOG,"Gwnxx Gwnyy Gwnzz Gwnxy Gwnxz Gwnyz ");				// Orientation tensors
		fprintf(TIMELOG,"Gwsxx Gwsyy Gwszz Gwsxy Gwsxz Gwsyz ");
		fprintf(TIMELOG,"Gnsxx Gnsyy Gnszz Gnsxy Gnsxz Gnsyz ");
		fprintf(TIMELOG,"trawn trJwn trRwn ");			//trimmed curvature,
		fprintf(TIMELOG,"wwndnw wwnsdnwn Jwnwwndnw "); 	//kinematic quantities,
		fprintf(TIMELOG,"Vw Aw Jw Xw "); 			//miknowski measures,
		fprintf(TIMELOG,"Vn An Jn Xn\n"); 			//miknowski measures,
	//	fprintf(TIMELOG,"Euler Kn Jn An\n"); 			//miknowski measures,
	}
}


// Destructor
SubPhase::~SubPhase()
{
	delete [] TempID;
    if ( TIMELOG!=NULL ) { fclose(TIMELOG); }

}

void SubPhase::SetParams(double rhoA, double rhoB, double tauA, double tauB, double force_x, double force_y, double force_z, double alpha)
{
	Fx = force_x;
	Fy = force_y;
	Fz = force_z;
	rho_n = rhoA;
	rho_w = rhoB;
	nu_n = (tauA-0.5)/3.f;
	nu_w = (tauB-0.5)/3.f;
	gamma_wn = 5.796*alpha;
}


